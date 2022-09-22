#include "cals.h"

#include <cmath>
#include <omp.h>

#include "multi_ktensor.h"
#include "utils/utils.h"
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
  CalsReport cp_cals(const Tensor &X, KtensorQueue &kt_queue, CalsParams &params)
  {
    // TODO investigate TTV library for fancier way to do this print.
    DEBUG(
        cout << "START method=Concurrent_ALS" << endl;
        cout << "Tensor size: [ ";
        for (const auto &mode: X.get_modes()) cout << mode << " ";
        cout << "]" << endl;
    )

    CalsReport rep;
    DEBUG(cout << "OMP_NUM_THREADS: " << omp_get_max_threads() << endl;)
    DEBUG(cout << "Threads set: " << get_threads() << endl;)

    rep.tensor_rank = X.get_rank();
    rep.n_modes = X.get_n_modes();
    rep.modes = X.get_modes();
    rep.X_norm = X.norm();

    rep.iter = 0;
    rep.max_iter = params.max_iterations;
    rep.n_threads = get_threads();
    rep.buffer_size = params.buffer_size;
    rep.n_ktensors = 0;
    rep.ktensor_rank_sum = 0;
    rep.tol = params.tol;
    rep.cuda = params.cuda;
    rep.update_method = params.update_method;
    rep.line_search = params.line_search;
    rep.line_search_interval = params.line_search_interval;
    rep.line_search_step = params.line_search_step;

#if WITH_TIME
    auto timer_size = std::ceil(20.0 * kt_queue.size() / rep.buffer_size) * rep.max_iter;
    rep.als_times = Matrix(ALS_TIMERS::LENGTH, timer_size);
    rep.mode_times = Matrix(MODE_TIMERS::LENGTH * rep.n_modes, timer_size);
    rep.mttkrp_times = Matrix(MTTKRP_TIMERS::LENGTH * rep.n_modes, timer_size);
    auto &als_times = rep.als_times;
    auto &mode_times = rep.mode_times;
    auto &mttkrp_times = rep.mttkrp_times;
    rep.flops_per_iteration.resize(timer_size, 0llu);
    rep.cols.resize(timer_size, 0);
#endif

    Timer total_time;
    Timer init_time;
    AlsTimers als_timers;
    ModeTimers mode_timers;

    using ALS_TIMERS = cals::AlsTimers::TIMERS;

    total_time.start();

    assert(rep.modes.size() >= 3);
    total_time.start();
    init_time.start();

    Matrix ten_workspace{};
    ls::LineSearchParams ls_params{};
    if (rep.line_search)
    {
      if (rep.modes.size() == 3)
      {
        ten_workspace = Matrix(rep.modes[0], rep.modes[1] * rep.modes[2]);
        ls_params.cuda = rep.cuda;
        ls_params.step = rep.line_search_step;
      } else
      {
        std::cerr << "Line search supported only for 3D tensors." << std::endl;
        abort();
      }
    }

    MultiKtensor mkt(rep.modes, params.buffer_size);
    mkt.set_cuda(rep.cuda);
    mkt.set_line_search(rep.line_search);

    auto &registry = mkt.get_registry();

    // Vector needed to hold the keys in the registry, so that multithreaded update & error can happen.
    // OpenMP (4.5) does not yet support parallalelizing non range-based loops.
    // Check OpenMP 5.0 (gcc 9) for possible improvement.
    vector<int> keys_v;
    keys_v.reserve(rep.buffer_size);

    // Create empty matrix to store last G for error calculation
    Matrix G_last(rep.modes[rep.n_modes - 1], params.buffer_size);

    // Calculate and Allocate Workspace for intermediate KRP and Twostep
    vector<Matrix> workspace;
    vector<dim_t> modes = rep.modes;  // Create a copy to sort it
    std::sort(modes.begin(), modes.end(), std::greater<>());
    vector<dim_t> sizes(std::ceil(static_cast<float>(modes.size() - 1) / 2.0),
                      0);  // sizes holds the nrows of each workspace matrix for KRP
    sizes[0] = modes[0] * modes[1];  // Two largest modes (Can be used as a workspace for twostep)
    for (auto i = 1lu; i < sizes.size(); i++) sizes[i] = modes[i + 1] * sizes[i - 1];

    workspace.resize(sizes.size());
    for (auto i = 0lu; i < sizes.size(); i++) workspace[i] = Matrix(sizes[i], params.buffer_size);

    // Vector needed to hold the keys of ktensors to be removed from the registry.
    // Can't traverse the map and remove at the same time safely.
    vector<int> remove_ids;
    remove_ids.reserve(100);

    // Create and initialize MTTKRP params
    mttkrp::MttkrpParams mttkrp_params;
    mttkrp_params.method = params.mttkrp_method;
    mttkrp_params.cuda = params.cuda;
    mttkrp_params.krp_params.cuda = params.cuda;
    mttkrp_params.lut = params.mttkrp_lut;

    init_time.stop();

#if CUDA_ENABLED
    auto stream = cuda::create_stream();
#endif

    if (params.cuda)
    {
#if CUDA_ENABLED
      if (X.get_cudata() == nullptr)
      {
        X.allocate_cudata(X.get_n_elements());
        X.send_to_device_async(stream);
      }

      for (auto &mf : mkt.get_factors())
        mf.allocate_cudata(mf.get_max_n_elements());

      for (auto &w : workspace)
        w.allocate_cudata(w.get_n_elements());
#else
      std::cerr << "Not compiled with CUDA support" << std::endl;
      exit(EXIT_FAILURE);
#endif
    }

    Timer lut_timer;
    lut_timer.start();
    if (mttkrp_params.lut.keys_v.empty())
    {
      std::cerr << "LUT empty" << endl;
      mttkrp_params.lut = mttkrp::read_lookup_table(rep.modes, get_threads(), params.cuda, false);
    }
    lut_timer.stop();
    //std::cout << "Lookup table creation took " << lut_timer.get_time() << " seconds." << std::endl;
    //std::cout << "Initialization took " << init_time.get_time() << " seconds." << std::endl;

    if (params.cuda)
    {
#if CUDA_ENABLED
      cudaStreamSynchronize(stream);
#else
      std::cerr << "Not compiled with CUDA support" << std::endl;
      exit(EXIT_FAILURE);
#endif
    }

    bool converged{false};
    mkt.set_iters(0);
    do
    {
      mkt.set_iters(mkt.get_iters() + 1);
      rep.iter = mkt.get_iters();
      // cout << "Iteration: " << rep.iter << endl;

      als_timers[ALS_TIMERS::ITERATION].start();

      while (!kt_queue.empty())
      {
        try
        {
          Ktensor &ktensor = kt_queue.front();
          mkt.add(ktensor);
          rep.n_ktensors++;
          rep.ktensor_rank_sum += ktensor.get_components();
          kt_queue.pop();
        }
        catch (BufferFull &e)
        {
          break;
        }
      }
      keys_v.clear();
      for (auto &[key, val] : registry)
        keys_v.emplace_back(key);

#pragma omp parallel for  // NOLINT(openmp-use-default-none)
      for (auto i = 0lu; i < keys_v.size(); i++)
      {
        auto &val = registry.at(keys_v[i]);
        if (rep.line_search && !((val.ktensor.get_iters() + 1) % rep.line_search_interval))
        {
          // TODO: create a Ktensor copy function to perform this copy with a single call.
          for (auto j = 0; j < rep.n_modes; j++)
            val.ls_ktensor.get_factor(j).copy(val.ktensor.get_factor(j));
        }
      }

#if WITH_TIME
      rep.cols[rep.iter - 1] = mkt.get_factor(0).get_cols();
      rep.flops_per_iteration[rep.iter - 1] = 0;
#endif

      // Loop over the modes
      for (auto n = 0; n < rep.n_modes; n++)
      {
        mode_timers[MODE_TIMERS::MTTKRP].start();

        auto &G = mttkrp::mttkrp(X, mkt, workspace, n, mttkrp_params);

        mode_timers[MODE_TIMERS::MTTKRP].stop();

        // Save the last G matrix for error calculation
        if (n == rep.n_modes - 1)
        {
          als_timers[ALS_TIMERS::G_COPY].start();
          G_last.resize(G.get_rows(), G.get_cols()).copy(G);
          als_timers[ALS_TIMERS::G_COPY].stop();
        }

        mode_timers[MODE_TIMERS::UPDATE].start();

#pragma omp parallel for  // NOLINT(openmp-use-default-none)
        for (auto i = 0lu; i < keys_v.size(); i++)
        {
          auto &val = registry.at(keys_v[i]);
          ops::hadamard_but_one(val.gramians, n);

          if (params.update_method == update::UPDATE_METHOD::UNCONSTRAINED)
            update::update_factor_unconstrained(val.ktensor.get_factor(n), val.gramians[n]);
          else
            update::update_factor_non_negative_constrained(val.ktensor.get_factor(n), val.gramians[n],
                                                           val.ktensor.get_active_set(n));

          // Question: rep.iter or val.ktensor.get_iter?
          val.ktensor.normalize(n, rep.iter);

          ops::update_gramian(val.ktensor.get_factor(n), val.gramians[n]);
        }

        if (params.cuda)
        {
#if CUDA_ENABLED
          mkt.get_factor(n).send_to_device();
#else
          std::cerr << "Not compiled with CUDA support" << std::endl;
          exit(EXIT_FAILURE);
#endif
        }

        mode_timers[MODE_TIMERS::UPDATE].stop();

#if WITH_TIME
        rep.flops_per_iteration[rep.iter - 1] += mttkrp_params.flops;
        for (auto i = 0; i < MODE_TIMERS::LENGTH; i++)
          mode_times(n * MODE_TIMERS::LENGTH + i, rep.iter - 1) = mode_timers[i].get_time();
        for (auto i = 0; i < MTTKRP_TIMERS::LENGTH; i++)
          mttkrp_times(n * MTTKRP_TIMERS::LENGTH + i, rep.iter - 1) = mttkrp_params.mttkrp_timers[i].get_time();
#endif
      } // for mode

      // Compute the approximation error
      als_timers[ALS_TIMERS::ERROR].start();

#pragma omp parallel for  // NOLINT(openmp-use-default-none)
      for (auto i = 0lu; i < keys_v.size(); i++)
      {
        auto &val = registry.at(keys_v[i]);

        const auto &ktensor_G_last = Matrix(val.ktensor.get_last_factor().get_rows(), val.ktensor.get_components(),
                                            G_last.get_data() + (val.col - mkt.get_start()) * G_last.get_col_stride());
        cals::ops::hadamard_all(val.gramians);
        const auto error = error::compute_fast_error(rep.X_norm, val.ktensor.get_lambda(),
                                                     val.ktensor.get_last_factor(),
                                                     ktensor_G_last, val.gramians[0]);
        assert((val.ktensor.get_iters() == 1) ||
               (val.ktensor.get_approximation_error() - error >= -1e-5));  // Make sure error decreases
        val.ktensor.set_approximation_error(error);
        val.ktensor.calculate_new_fit(rep.X_norm);
      }

      for (auto i = 0lu; i < keys_v.size(); i++)
      {
        auto &val = registry.at(keys_v[i]);
        if (rep.line_search && !(val.ktensor.get_iters() % (rep.line_search_interval)))
          ls::line_search(val.ktensor, val.ls_ktensor, workspace[0], ten_workspace, val.gramians, X, rep.X_norm,
                          ls_params);
      }

      als_timers[ALS_TIMERS::ERROR].stop();

      remove_ids.clear();
      for (auto &[key, val] : registry)
        if (!params.always_evict_first)
        {
          if (!params.force_max_iter)
          {
            if (val.ktensor.get_fit_diff() < rep.tol || val.ktensor.get_iters() >= rep.max_iter)
              remove_ids.push_back(key);
            else
              val.ktensor.set_iters(val.ktensor.get_iters() + 1);
          } else if (val.ktensor.get_iters() >= rep.max_iter)
            remove_ids.push_back(key);
          else
            val.ktensor.set_iters(val.ktensor.get_iters() + 1);
        } else  // Evict leftmost model for experiments
        {
          int id = mkt.get_leftmost_id();
          if (id != -1)
            remove_ids.push_back(id);
          break;
        }

      for (auto &key : remove_ids)
        mkt.remove(key);

      // Compression
      als_timers.timers[ALS_TIMERS::DEFRAGMENTATION].start();
      mkt.compress();
      als_timers.timers[ALS_TIMERS::DEFRAGMENTATION].stop();

      als_timers[ALS_TIMERS::ITERATION].stop();

#if WITH_TIME
      for (auto i = 0; i < ALS_TIMERS::LENGTH; i++)
        als_times(i, rep.iter - 1) = als_timers[i].get_time();
#endif
      DEBUG(
          cout << "CONVERGENCE " << rep.iter << endl;
          for (auto &[key, val]: registry) cout << val.ktensor.get_approximation_error() << " ";
          cout << endl;
      )

      if (kt_queue.empty() && registry.empty()) converged = true;

    } while (!converged);

    // Remove any model that has not converged (and copy it back to the original ktensor)
    if (!registry.empty())
    {
      remove_ids.clear();
      for (auto &[key, val] : registry)
        remove_ids.push_back(key);
      for (auto &key : remove_ids)
        mkt.remove(key);
    }

#if CUDA_ENABLED
    cuda::destroy_stream(stream);
#endif

    total_time.stop();
    rep.total_time = total_time.get_time();
    // cout << "Computation time: " << total_time.get_time() << endl;

    DEBUG(cout << "done." << endl;)

    return rep;
  }
}