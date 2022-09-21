#include <iomanip>
#include "als.h"

#include "utils/error.h"
#include "utils/line_search.h"

#if CUDA_ENABLED

#include "cuda_utils.h"

#endif

using std::cout;
using std::endl;
using ALS_TIMERS = cals::AlsTimers::TIMERS;
using MODE_TIMERS = cals::ModeTimers::TIMERS;
using MTTKRP_TIMERS = cals::MttkrpTimers::TIMERS;

namespace cals
{
  AlsReport cp_als(const Tensor &X, Ktensor &ktensor, AlsParams &params)
  {
    // TODO investigate TTV library for fancier way to do this print.
    DEBUG(
        cout << "START method=Regular_ALS" << endl;
        cout << "Tensor size: [ ";
        for (auto mode : X.get_modes()) cout << mode << " ";
        cout << "]" << endl;
        cout << "Ktensor rank= " << ktensor.get_components() << endl;
    )

    AlsReport rep;

    rep.tensor_rank = X.get_rank();
    rep.n_modes = X.get_n_modes();
    rep.modes = X.get_modes();
    rep.X_norm = X.norm();

    rep.iter = 0;
    rep.max_iter = params.max_iterations;
    rep.n_threads = get_threads();
    rep.ktensor_id = ktensor.get_id();
    rep.ktensor_rank = ktensor.get_components();
    rep.tol = params.tol;
    rep.cuda = params.cuda;
    rep.update_method = params.update_method;
    rep.line_search = params.line_search;
    rep.line_search_interval = params.line_search_interval;
    rep.line_search_step = params.line_search_step;

#if WITH_TIME
    rep.als_times = Matrix(ALS_TIMERS::LENGTH, rep.max_iter);
    rep.mode_times = Matrix(MODE_TIMERS::LENGTH * rep.n_modes, rep.max_iter);
    rep.mttkrp_times = Matrix(MTTKRP_TIMERS::LENGTH * rep.n_modes, rep.max_iter);
    auto &als_times = rep.als_times;
    auto &mode_times = rep.mode_times;
    auto &mttkrp_times = rep.mttkrp_times;
#endif

    Timer total_time;
    AlsTimers als_timers;
    ModeTimers mode_timers;

    using ALS_TIMERS = cals::AlsTimers::TIMERS;

    total_time.start();

    assert(rep.modes.size() >= 3);

    Ktensor ls_ktensor{};
    Matrix ten_workspace{};
    ls::LineSearchParams ls_params{};
    if (rep.line_search)
    {
      if (rep.modes.size() == 3)
      {
        ls_ktensor = Ktensor(rep.ktensor_rank, rep.modes);
        ten_workspace = Matrix(rep.modes[0], rep.modes[1] * rep.modes[2]);
        ls_params.cuda = rep.cuda;
        ls_params.step = rep.line_search_step;
      } else
      {
        std::cerr << "Line search supported only for 3D tensors." << std::endl;
        abort();
      }
    }

    // Allocate and Initialize Gramians vector
    auto gramians = vector<Matrix>(rep.n_modes);
    for (auto n = 0; n < rep.n_modes; n++) gramians[n] = Matrix(ktensor.get_components(), ktensor.get_components());
    for (auto n = 0; n < rep.n_modes; n++) cals::ops::update_gramian(ktensor.get_factor(n), gramians[n]);

    // Create empty matrix to store last G for error calculation
    Matrix G_last(rep.modes[rep.n_modes - 1], ktensor.get_components());

    // Calculate and Allocate Workspace for intermediate KRP and Twostep
    vector<Matrix> workspace;
    vector<dim_t> modes = rep.modes;  // Create a copy to sort it
    std::sort(modes.begin(), modes.end(), std::greater<>());
    vector<dim_t> sizes(std::ceil(static_cast<float>(modes.size() - 1) / 2.0),
                        0);  // sizes holds the nrows of each workspace matrix for KRP
    sizes[0] = modes[0] * modes[1];  // Two largest modes (Can be used as a workspace for twostep)
    for (auto i = 1lu; i < sizes.size(); i++)
      sizes[i] = modes[i + 1] * sizes[i - 1];

    workspace.resize(sizes.size());
    for (auto i = 0lu; i < sizes.size(); i++)
      workspace[i] = Matrix(sizes[i], ktensor.get_components());

    // Create and initialize MTTKRP params
    cals::mttkrp::MttkrpParams mttkrp_params;
    mttkrp_params.method = params.mttkrp_method;
    mttkrp_params.lut = params.mttkrp_lut;
    mttkrp_params.cuda = params.cuda;
    mttkrp_params.krp_params.cuda = params.cuda;

    Timer lut_timer;
    lut_timer.start();
    if (mttkrp_params.lut.keys_v.empty())
    {
      if (!params.suppress_lut_warning)
        std::cerr << "LUT empty" << endl;
      mttkrp_params.lut = mttkrp::read_lookup_table(rep.modes, get_threads(), params.cuda, params.suppress_lut_warning);
    }

    lut_timer.stop();
//    std::cout << "Lookup table creation took " << lut_timer.get_time() << "seconds." << std::endl;

#if CUDA_ENABLED
    auto stream = cuda::create_stream();
#endif

    if (params.cuda)
    {
#if CUDA_ENABLED
      if (!params.omp_enabled)
      {
        if (X.get_cudata() == nullptr)
        {
          X.allocate_cudata(X.get_n_elements());
          X.send_to_device();
        }
      }

      for (int i = 0; i < ktensor.get_n_modes(); i++)
      {
        ktensor.get_factor(i).allocate_cudata(ktensor.get_factor(i).get_n_elements());
        ktensor.get_factor(i).send_to_device();
      }

      for (auto &w : workspace)
        w.allocate_cudata(w.get_n_elements());
#else
      std::cerr << "Not compiled with CUDA support" << std::endl;
      exit(EXIT_FAILURE);
#endif
    }

    bool converged{false};
    rep.iter = 0;
    do
    {
      rep.iter++;
      rep.flops_per_iteration = 0;

      als_timers[ALS_TIMERS::ITERATION].start();

      if (rep.line_search && !((rep.iter + 1) % rep.line_search_interval))
      {
        for (auto i = 0; i < rep.n_modes; i++)
          ls_ktensor.get_factor(i).copy(ktensor.get_factor(i));
      }

      // Loop over the modes
      for (auto n = 0; n < rep.n_modes; n++)
      {
        auto &factor = ktensor.get_factor(n);

        mode_timers[MODE_TIMERS::MTTKRP].start();

        auto &G = mttkrp::mttkrp(X, ktensor, workspace, n, mttkrp_params);

        mode_timers[MODE_TIMERS::MTTKRP].stop();

        auto &H = cals::ops::hadamard_but_one(gramians, n);

        // Save the last G matrix for error calculation
        if (n == rep.n_modes - 1)
        {
          als_timers[ALS_TIMERS::G_COPY].start();
          G_last.copy(G);
          als_timers[ALS_TIMERS::G_COPY].stop();
        }

        mode_timers[MODE_TIMERS::UPDATE].start();

        if (params.update_method == cals::update::UPDATE_METHOD::UNCONSTRAINED)
          cals::update::update_factor_unconstrained(factor, H);
        else
          cals::update::update_factor_non_negative_constrained(factor, H, ktensor.get_active_set(n));

        ktensor.normalize(n, rep.iter);

        if (params.cuda)
        {
#if CUDA_ENABLED
          ktensor.get_factor(n).send_to_device_async(stream);
#else
          std::cerr << "Not compiled with CUDA support" << std::endl;
          exit(EXIT_FAILURE);
#endif
        }

        cals::ops::update_gramian(factor, gramians[n]);

        if (params.cuda)
        {
#if CUDA_ENABLED
          cudaStreamSynchronize(stream);
#else
          std::cerr << "Not compiled with CUDA support" << std::endl;
          exit(EXIT_FAILURE);
#endif
        }

        mode_timers[MODE_TIMERS::UPDATE].stop();

#if WITH_TIME
        rep.flops_per_iteration += mttkrp_params.flops;
        for (auto i = 0; i < MODE_TIMERS::LENGTH; i++)
          mode_times(n * MODE_TIMERS::LENGTH + i, rep.iter - 1) = mode_timers[i].get_time();
        for (auto i = 0; i < MTTKRP_TIMERS::LENGTH; i++)
          mttkrp_times(n * MTTKRP_TIMERS::LENGTH + i, rep.iter - 1) = mttkrp_params.mttkrp_timers[i].get_time();
#endif
      } // for mode

      // Compute the approximation error
      als_timers[ALS_TIMERS::ERROR].start();

      ktensor.set_iters(rep.iter);

      cals::ops::hadamard_all(gramians);
      const auto error = cals::error::compute_fast_error(rep.X_norm, ktensor.get_lambda(), ktensor.get_last_factor(),
                                                         G_last, gramians[0]);

      const auto old_error = ktensor.get_approximation_error();
      // Ensure error always decreases (first iteration excluded due to undefined initial error)
      const auto error_diff = old_error - error;
      if (rep.iter != 1 && (error_diff < -1e-8))
        std::cerr << std::scientific
                  << "error incr Ktensor: " << ktensor.get_id() << " Iter: " << std::setw(3) << rep.iter
                  << " || Error: " << std::setw(5) << error
                  << " || Old Error: " << std::setw(5) << old_error
                  << " || Diff: " << std::setw(5) << error_diff << endl;

      ktensor.set_approximation_error(error);
      ktensor.calculate_new_fit(rep.X_norm);

      als_timers[ALS_TIMERS::ERROR].stop();

      if (rep.line_search && !(rep.iter % rep.line_search_interval))
        ls::line_search(ktensor, ls_ktensor, workspace[0], ten_workspace, gramians, X, rep.X_norm, ls_params);

      als_timers[ALS_TIMERS::ITERATION].stop();

#if WITH_TIME
      for (auto i = 0; i < ALS_TIMERS::LENGTH; i++)
        als_times(i, rep.iter - 1) = als_timers[i].get_time();
#endif

      DEBUG(cout << "CONVERGENCE " << rep.iter << " " << ktensor.get_approximation_error() << endl;)

      if (!params.force_max_iter)
        converged = (ktensor.get_fit_diff() < rep.tol) || (rep.iter >= rep.max_iter);
      else
        converged = rep.iter >= rep.max_iter;

      if (params.force_max_iter) converged = (rep.iter >= rep.max_iter);
    } while (!converged);

    total_time.stop();
    rep.total_time = total_time.get_time();

#if CUDA_ENABLED
    cuda::destroy_stream(stream);
#endif

    DEBUG(cout << "done." << endl;)

    return rep;
  }

  std::vector<AlsReport> cp_omp_als(const Tensor &X, std::vector<Ktensor> &ktensor_v, AlsParams &params)
  {
    cals::Timer als_omp_time;
    als_omp_time.start();

    const auto n_ktensors = ktensor_v.size();

#if CUDA_ENABLED
    X.allocate_cudata(X.get_n_elements());
    X.send_to_device();
    params.omp_enabled = true;
#endif

    vector<AlsReport> reports(n_ktensors);
#pragma omp parallel for
    for (auto i = 0lu; i < n_ktensors; i++)
      reports[i] = cals::cp_als(X, ktensor_v[i], params);
    als_omp_time.stop();

#if CUDA_ENABLED
    params.omp_enabled = false;
#endif

    for (auto &rep : reports)
      rep.total_time = als_omp_time.get_time();

    return reports;
  }
}
