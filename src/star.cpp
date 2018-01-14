/*
 *  star.cpp
 *
 *  Created by Linar Mikeev (mikeev@cs.uni-saarland.de).
 *  Copyright (C) 2015 Saarland University. All rights reserved.
 *
 *
 *  This file is part of star.
 *
 *  star is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  star is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with star.  If not, see <http://www.gnu.org/licenses/>.
 */

#include <cstdlib>
#include <cstring>
#include <iostream>
#include <sstream>
#include <cfloat>
#include <cmath>
#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>
#include "experiment.hpp"
#include "../config.hpp"

#ifndef WIN32
#include <signal.h>
#include <sys/resource.h>
#endif

#if HAVE_MCR
#include "solver/matlab/libmatlabsd.hpp"
#endif

unsigned int experiment_id;

#ifdef STAR_WEB_INTERFACE

mysqlconnector dbconnector;

void sigabrt_handler(int) {
  if (!dbconnector.end_experiment(experiment_id, EXPERIMENT_STATUS_COMPLETED)) {
    std::cerr << "Error stopping experiment #" << experiment_id << std::endl;
  }
}

void sigterm_handler(int) {
  if (!dbconnector.end_experiment(experiment_id, EXPERIMENT_STATUS_STOPPED)) {
    std::cerr << "Error stopping experiment #" << experiment_id << std::endl;
  }
  exit(EXIT_SUCCESS);
}

void sigsegv_handler(int) {
  const char* err = "SIGSEGV";
  if (!dbconnector.end_experiment(experiment_id, EXPERIMENT_STATUS_FAILED,
                                  err)) {
    std::cerr << "Error stopping experiment #" << experiment_id << std::endl;
  }
  std::cerr << err << std::endl;
  exit(EXIT_FAILURE);
}

void sigxcpu_handler(int) {
  const char* err = "time limit exceeded";
  if (!dbconnector.end_experiment(experiment_id, EXPERIMENT_STATUS_FAILED,
                                  err)) {
    std::cerr << "Error stopping experiment #" << experiment_id << std::endl;
  }
  std::cerr << err << std::endl;
  exit(EXIT_FAILURE);
}

void sigxfsz_handler(int) {
  const char* err = "dump size limit exceeded";
  if (!dbconnector.end_experiment(experiment_id, EXPERIMENT_STATUS_FAILED,
                                  err)) {
    std::cerr << "Error stopping experiment #" << experiment_id << std::endl;
  }
  std::cerr << err << std::endl;
  exit(EXIT_FAILURE);
}

#endif

int run_main(int argc, char const** argv) {
  namespace pt = boost::property_tree;
  pt::ptree propTree;

  if (!load_config(propTree)) {
    return EXIT_FAILURE;
  }

#ifdef STAR_WEB_INTERFACE
  unsigned int user_id = 0;
  bool stop = false;
#endif

  namespace po = boost::program_options;
  po::options_description desc("Allowed options");
  desc.add_options()("help,h", "print usage message")
#ifdef STAR_WEB_INTERFACE
      ("experiment_id,e", po::value<unsigned int>(), "experiment id")(
          "user_id,u", po::value<unsigned int>(), "user id")(
          "stop", po::bool_switch(&stop), "stop experiment")
#else
      ("cleanout", "clean output")("experiment-file,e",
                                   po::value<std::string>(),
                                   "experiment file name")(
          "experiment-source,s", po::value<std::string>(), "experiment source")
#endif
      ;

  po::variables_map vm;
  try {
    store(parse_command_line(argc, argv, desc), vm);
  } catch (std::exception& e) {
    std::cerr << "error: " << e.what() << "\n";
    return EXIT_FAILURE;
  }

  if (vm.count("help")) {
    std::cout << desc << "\n";
    return EXIT_SUCCESS;
  }

#if HAVE_MCR
  char const* args[] = {"-nodesktop", "-nodisplay", "-nosplash", "-nojvm"};
  if (!mclInitializeApplication(args, 4)) {
    fprintf(stderr, "mclInitializeApplication = false\n");
    return EXIT_FAILURE;
  }
  if (!libmatlabsdInitialize()) {
    fprintf(stderr, "libmatlabsdInitialize = false\n");
    return EXIT_FAILURE;
  }
#endif

#ifdef STAR_WEB_INTERFACE

  if (!vm.count("experiment_id") || !vm.count("user_id")) {
    std::cerr << "Please specify experiment_id and user_id" << std::endl;
    return EXIT_FAILURE;
  }

  experiment_id = vm["experiment_id"].as<unsigned int>();
  user_id = vm["user_id"].as<unsigned int>();

  if (!db_connect(dbconnector, propTree)) {
    return EXIT_FAILURE;
  }

  if (vm.count("stop") && vm["stop"].as<bool>()) {
    if (!dbconnector.stop_experiment(experiment_id)) {
      return EXIT_FAILURE;
    }
    return EXIT_SUCCESS;
  }

  if (!dbconnector.is_admin(user_id)) {
#ifndef WIN32
    rlimit64 as_limit;
    if (getrlimit64(RLIMIT_AS, &as_limit) == 0) {
      as_limit.rlim_cur = 0x040000000LL * propTree.get("LIMITS.as", 64);
      if (setrlimit64(RLIMIT_AS, &as_limit) != 0) {
        std::cerr << strerror(errno) << std::endl;
      }
    }

    rlimit64 cpu_limit;
    if (getrlimit64(RLIMIT_CPU, &cpu_limit) == 0) {
      cpu_limit.rlim_cur = 60 * 60 * propTree.get("LIMITS.cpu", 336);
      if (setrlimit64(RLIMIT_CPU, &cpu_limit) != 0) {
        std::cerr << strerror(errno) << std::endl;
      }
    }

    rlimit64 fsize_limit;
    if (getrlimit64(RLIMIT_FSIZE, &fsize_limit) == 0) {
      fsize_limit.rlim_cur = 0x040000000LL * propTree.get("LIMITS.fsize", 1);
      if (setrlimit64(RLIMIT_FSIZE, &fsize_limit) != 0) {
        std::cerr << strerror(errno) << std::endl;
      }
    }
#endif
  }

  signal(SIGABRT, sigabrt_handler);
  signal(SIGTERM, sigterm_handler);
  signal(SIGSEGV, sigsegv_handler);
  signal(SIGXCPU, sigxcpu_handler);
  signal(SIGXFSZ, sigxfsz_handler);

  experiment E;
  if (!E.run(&dbconnector, user_id, experiment_id)) {
    return EXIT_FAILURE;
  }

#else

  if (vm.count("cleanout")) {
    namespace fs = boost::filesystem;

    fs::path output_path("output");
    for (fs::directory_iterator end_dir_it, it(output_path); it != end_dir_it;
         ++it) {
      fs::remove_all(it->path());
    }

    return EXIT_SUCCESS;
  }

  std::string experiment_file;
  if (vm.count("experiment-file")) {
    experiment_file = vm["experiment-file"].as<std::string>();
  }

  std::string experiment_source;
  if (vm.count("experiment-source")) {
    experiment_source = vm["experiment-source"].as<std::string>();
  }

  experiment E;
  if (experiment_file != "") {
    if (experiment_source != "") {
      std::cerr << "Please specify either experiment file or experiment source"
                << std::endl;
      return EXIT_FAILURE;
    }

    if (!E.run_file(experiment_file)) {
      return EXIT_FAILURE;
    }
  } else if (experiment_source != "") {
    if (experiment_file != "") {
      std::cerr << "Please specify either experiment file or experiment source"
                << std::endl;
      return EXIT_FAILURE;
    }

    if (!E.run_src(experiment_source)) {
      return EXIT_FAILURE;
    }
  } else {
    std::cerr << "Please specify either experiment file or experiment source"
              << std::endl;
    return EXIT_FAILURE;
  }
#endif

  return EXIT_SUCCESS;
}

int main(int argc, char const** argv) {
#ifdef HAVE_MCR
  mclmcrInitialize();
  return mclRunMain((mclMainFcnType)run_main, argc, argv);
#else
  return run_main(argc, argv);
#endif
}

#ifndef STAR_WEB_INTERFACE
int star_run_file(lua_State*, const char* experiment_fname, star_fnsi* fnsi) {
  experiment E(fnsi);
  return E.run_file(experiment_fname);
}

int star_run_src(lua_State*, const char* experiment_src, star_fnsi* fnsi) {
  experiment E(fnsi);
  return E.run_src(experiment_src);
}
#endif
