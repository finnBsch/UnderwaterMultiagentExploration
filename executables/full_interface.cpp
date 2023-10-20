//
// Created by finn on 7/2/23.
//
#include <sys/stat.h>
#include "sim_interface.h"
int main(int argc, char **argv){
    std::vector<std::string> names = {};
    std::vector<int> repeats = {};
    std::vector<int> num_agents = {};

    if (argc > 1) {
        names.push_back(argv[1]);
        num_agents.push_back(std::stoi(argv[2]));
        repeats.push_back(std::stoi(argv[3]));
    }
    else{
        names = {"CaseA", "CaseB", "CaseC", "CaseA2", "CaseA1"};
        repeats = {0};
        num_agents = {4, 3, 2, 1};
    }
    std::string out_path = "../../"
                           "/Experiments/ClosedLoop/NumberOfAgents/TimingTest/";

    for(int i = 0; i < names.size(); i++) {
        for (auto j: num_agents) {
            std::string log_base_path =
                    out_path + names[i] + "/" + std::to_string(j) + "/";
            struct stat sb;
            if (stat(log_base_path.c_str(), &sb) == 0) {
            } else {
                mkdir(log_base_path.c_str(), 0777);
                std::cout << "Created Experiment Directory in "
                          << log_base_path << std::endl;
            }
        }

        for (auto j : num_agents) {
            for (auto k : repeats) {
                ConfigFileScenario scenario(names[i]);
                GmrfParams gmrf_params;
                RRTParams rrt_params;
                SQPParams sqp_params;
                sqp_params.num_agents = j ;
                FullInterface full(out_path + names[i] + "/" + std::to_string(j ) + "/" + std::to_string(k) + "/",
                                   &sqp_params,
                                   &rrt_params,
                                   &gmrf_params, &scenario,
                                   StoppingCritereon::N_steps, true, false,
                                   40.0,
                                   2.0);
                try {
                    full.runFull();
                }
                catch (std::exception &e) {
                    std::cout << "Exception: " << e.what() << std::endl;
                }
            }
        }
    }
    return 0;
}