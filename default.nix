let
  pkgs_stable = import <nixos> {};
  pkgs = import <nixos-unstable> {};
  crlibm = import ./crlibm.nix {
    pkgs = pkgs;
    stdenv = pkgs.gcc10Stdenv;
  };
in
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
        # profiling and debugging
        gdb
        valgrind
        linuxPackages.perf
        cmake
        eigen
        doctest
        clang_10
        openlibm
        crlibm
        mpfi
        pkgconfig
        pkgs_stable.julia_11
        # visualize profile results
        # qcachegrind
      ];
    }
