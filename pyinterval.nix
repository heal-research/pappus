{ pkgs, lib, stdenv, buildPythonPackage, fetchPypi }:
let
  crlibm = import ./crlibm.nix {
    pkgs = pkgs;
    stdenv = pkgs.gcc10Stdenv;
  };
  pycrlibm = import ./pycrlibm.nix { inherit pkgs lib stdenv buildPythonPackage fetchPypi; };
in
  buildPythonPackage rec {
    pname = "pyinterval";
    version = "1.2.0";

    buildInputs = [ (with pkgs.python38Packages; [ six pycrlibm ]) ];

    src = fetchPypi {
      inherit pname version;
      sha256 = "8c46224a05815affa803ed5620432236bba620a246701b860d98ae7358a41d21";
    };

    doCheck = false;

    meta = with lib; {
      homepage = https://github.com/taschini/pyinterval;
      description = "Python implementation of an algebraically closed interval system on the extended real number set.";
    };
  }
