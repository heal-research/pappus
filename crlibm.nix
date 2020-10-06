{ pkgs, stdenv }:

stdenv.mkDerivation rec {
  pname = "crlibm";
  version = "1.0beta5";

  nativeBuildInputs = with pkgs; [ gcc10 autoreconfHook pkgconfig ];

  CFLAGS = "-march=znver2 -O3";

  src = pkgs.fetchgit {
    url = "https://scm.gforge.inria.fr/anonscm/git/metalibm/crlibm.git";
    rev = "eb3063791aa75bc9705b49283bf14250465220a7";
    sha256 = "0dlwa26vji9c3cwj5nc58dh857pa22x00cpn7l3rajnvarar6f0x";
  };

  installPhase = ''
      ./prepare --prefix=/
      make
      make DESTDIR=$out install
  '';

  meta = with stdenv.lib; {
    description = "An efficient and proven correctly-rounded mathematical library.";
  };
}


