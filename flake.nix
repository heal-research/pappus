{
  description = "pappus affine arithmetic library";

  inputs = {
    flake-parts.url = "github:hercules-ci/flake-parts";
    nixpkgs.url = "github:nixos/nixpkgs/master";
    foolnotion.url = "github:foolnotion/nur-pkg";
    foolnotion.inputs.nixpkgs.follows = "nixpkgs";
    aaflib.url = "github:foolnotion/aaflib";
    aaflib.inputs.nixpkgs.follows = "nixpkgs";
    aaflib.inputs.foolnotion.follows = "foolnotion";
    aaflib.inputs.flake-parts.follows = "flake-parts";
  };

  outputs =
    inputs@{
      self,
      flake-parts,
      foolnotion,
      nixpkgs,
      aaflib,
    }:
    flake-parts.lib.mkFlake { inherit inputs; } rec {
      systems = [
        "x86_64-linux"
        "x86_64-darwin"
        "aarch64-linux"
        "aarch64-darwin"
      ];

      perSystem =
        { pkgs, system, inputs', ... }:
        let
          pkgs = import nixpkgs {
            inherit system;
            overlays = [ foolnotion.overlay ];
          };
          stdenv_ = pkgs.llvmPackages_latest.stdenv;

          pappus = stdenv_.mkDerivation {
            name = "pappus";
            src = self;

            cmakeFlags = [
              "-DBUILD_TESTS=OFF"
              "-DUSE_COREMATH=OFF"
              "-DCMAKE_BUILD_TYPE=Release"
            ];

            nativeBuildInputs = with pkgs; [ cmake ];

            buildInputs = with pkgs; [
              gch-small-vector
              eve
              openlibm
            ];
          };

          enableTesting = true;
        in
        rec {
          packages.default = pappus;

          devShells.default = stdenv_.mkDerivation {
            name = "pappus-dev";

            nativeBuildInputs =
              pappus.nativeBuildInputs
              ++ (with pkgs; [
                clang-tools
                bear
              ]);

            buildInputs =
              pappus.buildInputs
              ++ (
                with pkgs;
                pkgs.lib.optionals pkgs.stdenv.isLinux [
                  gdb
                  valgrind
                  linuxPackages.perf
                ]
              )
              ++ pkgs.lib.optionals enableTesting [
                pkgs.catch2_3
                pkgs.nanobench
                inputs'.aaflib.packages.default
              ];
          };
        };
    };
}
