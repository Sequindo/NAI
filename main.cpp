#include <iostream>
#include <vector>
#include <math.h>
#include <cstdlib>
#include <functional>
#include <random>
#include <string>
#define RASTRIGRIN_MIN -5.121
#define RASTRIGRIN_MAX 5.121
#define HIMMELBLAU_MIN -5
#define HIMMELBLAU_MAX 5

#define DEFAULT_ITER_VALUE 10000

std::random_device rd;  //Will be used to obtain a seed for the random number engine
std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()

int main(int argc, char** argv)
{
    std::uniform_real_distribution<> distrib_r(HIMMELBLAU_MIN, HIMMELBLAU_MAX);
    std::vector<double> p0 = {
        distrib_r(gen),
        distrib_r(gen),
    };

    auto rastrigrin_domain = [](std::vector<double>args)
    {
        bool is_domain = true;
        for(auto &arg : args)
        {
            is_domain = is_domain&&(arg > RASTRIGRIN_MIN && arg < RASTRIGRIN_MAX);
        }
        return is_domain;
    };
    auto rastrigrin = [](std::vector<double>args)
    {
        const int A = 10;
        int n = args.size();

        double res = A*n;
        for(auto &arg : args)
        {
            res = res + arg*arg - A*cos(2*M_PI*arg);
        }
        return res;
    };

    auto himmelblau_domain = [](std::vector<double>args)
    {
        bool is_domain = true;
        for(auto arg : args)
        {
            is_domain = is_domain&&(arg >= HIMMELBLAU_MIN && arg <= HIMMELBLAU_MAX);
        }
        return is_domain;
    };
    auto himmelblau = [](std::vector<double>args)
    {
        double x = args.at(0); 
        double y = args.at(1);
        return pow((x*x+y-11),2) + pow((x+y*y-7),2);
    };

    auto generate_gnuplot_output = [](double border_min, double step, auto func) //generate gnuplot output for whole function domain
    {
        std::vector<double> args;
        double x, y = (fabs(border_min))*-1.0;;
        double border_max = fabs(border_min);
        while(y<border_max)
        {   
            x = (fabs(border_min))*-1.0;;
            while(x<border_max)
            {
                args.push_back(x);
                args.push_back(y);
                x = x + step;
            }
            y = y + step;
        }
        while(!args.empty())
        {
            auto y_arg = args.back();
            args.pop_back();
            auto x_arg = args.back();
            args.pop_back();
            std::cout << x_arg << " " << y_arg << " " << func(std::vector<double>{x_arg, y_arg}) << std::endl;
        }
    };
    //generate_gnuplot_output(RASTRIGRIN_MIN, 0.05, rastrigrin); 

    auto hill_climb = [](std::function<double(std::vector<double>)> function, std::function<bool(std::vector<double>)> f_domain,
                        std::vector<double> p0, int iterations, double delta)
    {
        if(!f_domain(p0))
        {
            throw std::string("Points "+std::to_string(p0.at(0))+", "+std::to_string(p0.at(1))+" are out of function's domain");
        }
        auto p = p0;
        std::uniform_int_distribution<> distrib(0, p.size() - 1); //choose vector element
        delta = fabs(delta);
        std::uniform_real_distribution<> distrib_r((-1.0*delta), delta);
        for (int i = 0; i < iterations; i++) {
            auto p2 = p;

            p2[distrib(gen)] += distrib_r(gen);
            if(!f_domain(p2)) continue;
            double y2 = function(p2);
            if (y2 < function(p)) {
                p = p2;
            }
        }
        return p;
    };

    auto simulated_annealing = [](std::function<double(std::vector<double>)> function, std::function<bool(std::vector<double>)> f_domain,
                        std::vector<double> p0, int iterations, std::function<std::vector<double>(std::vector<double>)> neighbours, std::function<double(int)> temp)
    {
        auto s_current = p0;
        auto s_global_best = p0;

        std::uniform_real_distribution<> u_k(0.0, 1.0);

        for(int i=0;i<iterations;i++)
        {
            auto s_next = neighbours(s_current);
            if(function(s_next) < function(s_current))
            {
                s_current = s_next;
            }
            else
            {
                double u = u_k(gen);
                if(u < exp(-abs(function(s_next) - function(s_current)) / temp(i)))
                {
                    s_current = s_next;
                }
            }
            if(s_current > s_global_best)
            {
                s_global_best = s_current;
            }
        }
    };
    std::function<bool(std::vector<double>)> domain_func = himmelblau_domain;
    std::function<double(std::vector<double>)> optimalization_func = himmelblau;
    auto iterations_num = DEFAULT_ITER_VALUE;
    if(argc>2)
    {
        if(argv[1]=="rastrigrin")
        {
            domain_func = rastrigrin_domain;
            optimalization_func = rastrigrin;
        }
        iterations_num = std::stoi(argv[2]);
    }
    
    double ss = 0.01;

    for(int it=0;it<1000;it++)
    {
        p0 = {distrib_r(gen), distrib_r(gen)};
        try
        {
            auto res = hill_climb(optimalization_func, domain_func, p0, iterations_num, ss);
            while(!res.empty())
            {
                auto y_arg = res.back();
                res.pop_back();
                auto x_arg = res.back();
                res.pop_back();
                std::cout << x_arg << " " << y_arg << " " << optimalization_func(std::vector<double>{x_arg, y_arg}) << std::endl;
            }
        }
        catch(const std::string& str)
        {
            //std::cerr << str << '\n';
        } 
    }
    return 0;
}