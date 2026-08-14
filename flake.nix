{
  description = "amplicon-qc - FASTQ read counter";

  inputs = {
    nixpkgs.url = "github:NixOS/nixpkgs/nixos-25.05";
    flake-utils.url = "github:numtide/flake-utils";
  };

  outputs = {
    self,
    nixpkgs,
    flake-utils,
  }:
    flake-utils.lib.eachDefaultSystem (system: let
      pkgs = import nixpkgs {inherit system;};

      version =
        if self ? rev
        then builtins.substring 0 8 self.rev
        else "dev";
    in {
      packages.default = pkgs.rustPlatform.buildRustPackage {
        pname = "amplicon-qc";
        inherit version;

        src = self;

        cargoLock = {
          lockFile = ./Cargo.lock;
        };

        nativeBuildInputs = with pkgs; [
          pkg-config
        ];

        meta = with pkgs.lib; {
          description = "FASTQ flank check tool";
          platforms = platforms.linux;
        };
      };

      defaultPackage = self.packages.${system}.default;

      apps.default = {
        type = "app";
        program = "${self.packages.${system}.default}/bin/amplicon-qc";
      };

      devShells.default = pkgs.mkShell {
        packages = with pkgs; [
          cargo
          rustc
          rustfmt
          clippy
          rust-analyzer
          pkg-config
        ];
      };
    });
}
