{
  inputs = {
    nixpkgs.url = "github:NixOS/nixpkgs/release-25.11"; # that's 23.05
    utils.url = "github:numtide/flake-utils";
    naersk.url = "github:nmattia/naersk";
    naersk.inputs.nixpkgs.follows = "nixpkgs";
    rust-overlay.url = "github:oxalica/rust-overlay";
    rust-overlay.inputs.nixpkgs.follows = "nixpkgs";
  };

  outputs =
    {
      self,
      nixpkgs,
      utils,
      naersk,
      rust-overlay,
    }:
    utils.lib.eachDefaultSystem (
      system:
      let
        #pkgs = nixpkgs.legacyPackages."${system}";
        overlays = [ (import rust-overlay) ];
        pkgs = import nixpkgs { inherit system overlays; };
        rust = pkgs.rust-bin.stable."1.92.0".default.override {
          targets = [ "x86_64-unknown-linux-musl" ];
          extensions = [ "llvm-tools-preview"  "rust-analyzer"];
        };

        # Override the version used in naersk
        naersk-lib = naersk.lib."${system}".override {
          cargo = rust;
          rustc = rust;
        };

        bacon = pkgs.bacon;
      in
      rec {
        # `nix build`
        packages.mbf-fastq-processor = naersk-lib.buildPackage {
          pname = "mbf-fastq-processor";
          root = ./.;
          nativeBuildInputs = with pkgs; [
            pkg-config
            cmake
            gcc
            gnumake
          ];
          buildInputs = with pkgs; [
            openssl
            rapidgzip
            which
            mold
          ];
          release = true;
          CARGO_PROFILE_RELEASE_debug = "0";
          COMMIT_HASH = self.rev or (pkgs.lib.removeSuffix "-dirty" self.dirtyRev or "unknown-not-in-git");
          NIX_RAPIDGZIP = "${pkgs.rapidgzip}/bin/rapidgzip";

          # copyBinsFilter = ''
          #   select(.reason == "compiler-artifact" and .executable != null and .profile.test == false and .target.name != "mbf-fastq-processor-test-runner")
          # '';
        };
        packages.mbf-fastq-processor_other_linux =
          (naersk-lib.buildPackage {
            pname = "mbf-fastq-processor";
            root = ./.;
            nativeBuildInputs = with pkgs; [
              pkg-config
              cmake
              gcc
              gnumake
              mold
            ];
            buildInputs = with pkgs; [ openssl ];
            release = true;
            CARGO_PROFILE_RELEASE_debug = "0";
            COMMIT_HASH = self.rev or (pkgs.lib.removeSuffix "-dirty" self.dirtyRev or "unknown-not-in-git");
            # copyBinsFilter = ''
            #   select(.reason == "compiler-artifact" and .executable != null and .profile.test == false and .target.name != "mbf-fastq-processor-test-runner")
            # '';
          }).overrideAttrs
            {
              # make it compatible with other linuxes. It's statically linked anyway
              postInstall = ''
                patchelf $out/bin/mbf-fastq-processor --set-interpreter "/lib64/ld-linux-x86-64.so.2"
              '';
            };
        packages.mbf-fastq-processor-docker =
          let
            binary = packages.mbf-fastq-processor_other_linux;
          in
          pkgs.dockerTools.buildLayeredImage {
            name = "mbf-fastq-processor";
            tag = "latest";
            # provide a minimal base with glibc and a busybox shell
            contents = [
              pkgs.busybox
              pkgs.glibc
              binary
            ];
            config = {
              Env = [ "PATH=/usr/local/bin:/bin" ];
              Entrypoint = [ "/bin/mbf-fastq-processor" ];
              WorkingDir = "/work";
            };
          };
        packages.check = naersk-lib.buildPackage {
          src = ./.;
          mode = "check";
          name = "mbf-fastq-processor";
          nativeBuildInputs = with pkgs; [
            pkg-config
            cmake
            gcc
            gnumake
            zstd.bin
            mold
          ];
          buildInputs = with pkgs; [ openssl ];
        };
        packages.test = naersk-lib.buildPackage {
          # not using naersk test mode, it eats the binaries, we need that binary
          pname = "mbf-fastq-processor";
          root = ./.;
          nativeBuildInputs = with pkgs; [
            pkg-config
            cmake
            gcc
            gnumake
            python3
            rapidgzip
            which
            mold
          ];
          buildInputs = with pkgs; [ openssl ];
          release = true;
          CARGO_PROFILE_RELEASE_debug = "0";
          COMMIT_HASH = self.rev or (pkgs.lib.removeSuffix "-dirty" self.dirtyRev or "unknown-not-in-git");
          NIX_RAPIDGZIP_ = "${pkgs.rapidgzip}/bin/rapidgzip"; # note the _, it's special cased.
          RUST_LOG = "trace";
          # every other test happens wit hteh rapidgzip in the path.
          postInstall = ''
            # run the friendly panic test, expect a non 0 return code.
            # capture stderr

            result=$( { cargo run --release --bin mbf-fastq-processor -- --test-friendly-panic 1>/dev/null; } 2>&1 ) || status=$? : "${"status:=0"}"
            if [ "$status" -eq 0 ]; then
              echo "Unexpected success"
              exit 1
            fi
            if [[ ! $result =~ "this is embarrassing" ]]; then
                echo "Error: friendly panic message ' not found in stderr"
                exit 1
            fi
            # without NIX_RAPIDGZIP, the test passes because the error is thrown
            echo 'without NIX_RAPIDGZIP'
            cargo test --release 
            # but with NIX_RAPIDGZIP, the test fails because there is a fallback

            echo 'with NIX_RAPIDGZIP'
            set +e  # Temporarily disable exit-on-error
            NIX_RAPIDGZIP=$NIX_RAPIDGZIP_ cargo test --release error_no_rapid_gzip
            set -e  # Re-enable exit-on-error
            if [ "$status" -eq 0 ]; then
              echo "Unexpected success when testing no-rapid-gzip-error-case"
              exit 1
            fi

          '';

          # src = ./.;
          # buildInputs = with pkgs; [openssl ];
          # mode = "test";
          # nativeBuildInputs = with pkgs; [pkg-config cargo-nextest];
          # cargoTestCommands = old: ["cargo nextest run $cargo_test_options --no-fail-fast"];
          # copySources = ["tests" "test_cases" "dev"];
          # copyBins = true;

          # override = {
          #   buildPhase = ":";
          #   postCheck = ''
          #      # make sure that the friendly panic test outputs a friendly panic
          #      ls -la
          #     cargo build --release
          #      if [ $? -ne 0 ]; then
          #          echo "Error: Command failed with non-zero status code"
          #          exit 1
          #      fi
          #      result=`cargo run --release -- --friendly-panic-test`

          #      # Check if stderr contains 'this is embarrasing'
          #      if grep -q "this is embarrasing" <(echo "$result"); then
          #          echo "Error: 'this is embarrasing' found in stderr"
          #          exit 1
          #      fi

          #      # now run our actual test cases
          #      cat Cargo.toml
          #     cargo run --release --bin mbf-fastq-processor-test-runner test_cases
          #   '';
          # };
          # doCheck = true;
        };
        # haven't been able to get this to work
        # packages.coverage = naersk-lib.buildPackage {
        #   src = ./.;
        #   buildInputs = with pkgs; [openssl ];
        #   mode = "test";
        #   nativeBuildInputs = with pkgs; [pkg-config cargo-nextest cargo-llvm-cov];
        #   cargoTestCommands = old: ["cargo llvm-cov nextest --no-tests=fail --run-ignored all"];
        #   override = {
        #     buildPhase = ":";
        #     postCheck = ''
        #       cp  target/llvm-cov/html $out/ -r
        #       '';
        #   };
        #   doCheck = true;
        # };
        #cargoTestCommands = old: ["cargo llvm-cov --html nextest --verbose $cargo_test_options"];

        defaultPackage = packages.mbf-fastq-processor;

        # `nix run`
        apps.mbf-fastq-processor = utils.lib.mkApp { drv = packages.mbf-fastq-processor; };
        defaultApp = apps.mbf-fastq-processor;

        # `nix develop`
        devShell = pkgs.mkShell {
          COMMIT_HASH = self.rev or (pkgs.lib.removeSuffix "-dirty" self.dirtyRev or "unknown-not-in-git");
          # we only link with mold in our dev environment for build speed. CI can use the old school rust linker
          shellHook = ''
            export RUSTFLAGS="-C link-arg=-fuse-ld=mold"
            # Set shell for cmake builds
            export CONFIG_SHELL="${pkgs.bash}/bin/bash"
            export SHELL="${pkgs.bash}/bin/bash"
          '';
          # supply the specific rust version
          nativeBuildInputs = [
            bacon
            pkgs.bash
            pkgs.cargo-audit
            pkgs.cargo-bloat
            pkgs.cargo-crev
            pkgs.cargo-deny
            pkgs.cargo-features-manager
            pkgs.cargo-flamegraph
            pkgs.cargo-insta
            pkgs.cargo-license
            pkgs.cargo-llvm-cov
            pkgs.cargo-llvm-lines
            pkgs.cargo-machete
            pkgs.cargo-mutants
            pkgs.cargo-nextest
            pkgs.cargo-outdated
            pkgs.cargo-shear
            #pkgs.cargo-udeps
            pkgs.cargo-vet
            pkgs.cmake
            pkgs.gcc
            pkgs.gnumake
            pkgs.git
            pkgs.hugo
            pkgs.jq
            pkgs.mold
            pkgs.openssl
            pkgs.pkg-config
            pkgs.samply
            (pkgs.python3.withPackages (
              ps: with ps; [
                scipy
                pysam
                toml
              ]
            ))
            pkgs.rapidgzip
            pkgs.which
            pkgs.ripgrep
            #rust.rust-analyzer
            pkgs.shellcheck
            rust
          ];
        };
      }
    );
}
# {
