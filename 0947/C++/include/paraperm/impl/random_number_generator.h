#ifndef PARPAPERM_IMPL_RANDOM_NUMBER_GENERATOR_H
#define PARPAPERM_IMPL_RANDOM_NUMBER_GENERATOR_H

#include <ctime>
#include <cstdint>
#include <cstdlib>

#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_int_distribution.hpp>
#include <boost/shared_ptr.hpp>

namespace paraperm
{
    namespace utils
    {
        template <typename T>
        struct random_number_generator_traits 
        {
            typedef boost::random::mt19937 engine_type;
            typedef boost::random::mt19937::result_type engine_result_type;
            typedef boost::random::uniform_int_distribution<T> dist_type;
        };

        template<>
        struct random_number_generator_traits<uint64_t>
        {
            typedef boost::random::mt19937_64 engine_type;
            typedef boost::random::mt19937::result_type engine_result_type;
            typedef boost::random::uniform_int_distribution<uint64_t> dist_type;
        };

        template <typename INT>
        class random_number_generator
        {
            public:
                typedef INT result_type;

       explicit random_number_generator(INT rank_)
                {
                    construct(0, std::numeric_limits<INT>::max(), rank_);
                }

                random_number_generator(INT a, INT b, INT rank_)
                {
                    construct(a, b, rank_);
                }

                INT operator()()
                {
                    return (*dist)(*rng);
                }

            private:
                template <typename T>
                struct traits;

                typedef typename random_number_generator_traits<INT>::engine_type engine_type;
                typedef typename random_number_generator_traits<INT>::engine_result_type engine_result_type;
                typedef typename random_number_generator_traits<INT>::dist_type dist_type;

                boost::shared_ptr<engine_type> rng;
                boost::shared_ptr<dist_type> dist;

                void construct(INT a, INT b, INT rank_)
                {
                    std::time_t seconds = std::time(NULL);
                    long temp = (seconds * 181) * ((rank_ - 83) * 359);
                    engine_result_type seed = std::abs(temp) % 104729;

                    rng.reset(new engine_type(seed));
                    dist.reset(new dist_type(a, b));
                }

        };

    }

}

#endif
