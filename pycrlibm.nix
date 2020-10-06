{ pkgs, lib, stdenv, buildPythonPackage, fetchPypi }:
let
  crlibm = import ./crlibm.nix { inherit pkgs stdenv; };
in
  buildPythonPackage rec {
    pname = "pycrlibm";
    version = "1.0.3";

    buildInputs = [ (with pkgs.python38Packages; [ six ]) crlibm ];

    src = fetchPypi {
      pname = "crlibm";
      inherit version;
      sha256 = "48e17981f90d69c6bb0013f68bacbe7a157de864a533d15dd196ca7e98348a35";
    };

    doCheck = false;

    meta = with lib; {
      homepage = https://github.com/taschini/pycrlibm;
      description = "Stupid wrapper.";
    };
  }
