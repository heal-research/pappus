let
  pkgs_stable = import <nixos> {};
  pkgs = import <nixos-unstable> {};
in
#unstable.llvmPackages_10.stdenv.mkDerivation {
pkgs.gcc10Stdenv.mkDerivation {
    name = "pappus-env";
    hardeningDisable = [ "all" ]; 

    buildInputs = with pkgs; [
        # python environment for bindings and scripting
        #(pkgs.python38.withPackages (ps: with ps; [ pip numpy pandas pybind11 pyperf colorama coloredlogs seaborn ]))
        (pkgs_stable.python38.withPackages (ps: with ps; [ sphinx recommonmark sphinx_rtd_theme ]))
        # Project dependencies
        ccls # completion vim
        bear # generate compilation database
        tbb
        # profiling and debugging
        gdb
        valgrind
        linuxPackages.perf
        git
        cmake
        eigen
        doctest
        clang_10
        # visualize profile results
        qcachegrind
      ];
    }
