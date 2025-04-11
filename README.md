# Condor-rs - post-quantum proof schemes and signatures

## Crates

* [Labrador](https://nethermindeth.github.io/condor-rs/labrador/) post-quantum proof system


## Development

Enable [commit signing](https://docs.github.com/en/authentication/managing-commit-signature-verification/signing-commits)

```sh
git config commit.gpgsign true
```

### Prerequisites

* [Rust](https://www.rust-lang.org/tools/install)
* [cargo deny](https://github.com/EmbarkStudios/cargo-deny)
* [typos](https://github.com/crate-ci/typos?tab=readme-ov-file#install)

### Code quality assurance

Install a pre-push git hook:

```sh
git config core.hooksPath .githooks
```

### Running the Rust Documentation Locally
After cloning the repository, follow the instructions below to run the documentation locally:

```sh
cargo doc
```

Docs for `labrador`:

```sh
RUSTDOCFLAGS="--html-in-header katex-header.html" cargo doc --no-deps -p labrador --open
```

## License

Apache 2.0

## Contributors ✨

Thanks goes to these wonderful people ([emoji key](https://allcontributors.org/docs/en/emoji-key)):

<!-- ALL-CONTRIBUTORS-LIST:START - Do not remove or modify this section -->
<!-- prettier-ignore-start -->
<!-- markdownlint-disable -->
<table>
  <tbody>
    <tr>
      <td align="center" valign="top" width="14.28%"><a href="https://github.com/frozenspider"><img src="https://avatars.githubusercontent.com/u/2077017?v=4?s=100" width="100px;" alt="Alex Abdugafarov"/><br /><sub><b>Alex Abdugafarov</b></sub></a><br /><a href="https://github.com/NethermindEth/condor-rs/commits?author=frozenspider" title="Code">💻</a></td>
      <td align="center" valign="top" width="14.28%"><a href="https://github.com/mattsuffern"><img src="https://avatars.githubusercontent.com/u/135047609?v=4?s=100" width="100px;" alt="Mateo Suffern"/><br /><sub><b>Mateo Suffern</b></sub></a><br /><a href="https://github.com/NethermindEth/condor-rs/commits?author=mattsuffern" title="Code">💻</a></td>
      <td align="center" valign="top" width="14.28%"><a href="https://github.com/pycckuu"><img src="https://avatars.githubusercontent.com/u/1489583?v=4?s=100" width="100px;" alt="Igor Markelov"/><br /><sub><b>Igor Markelov</b></sub></a><br /><a href="https://github.com/NethermindEth/condor-rs/commits?author=pycckuu" title="Code">💻</a></td>
      <td align="center" valign="top" width="14.28%"><a href="https://github.com/Yue-Zhou1"><img src="https://avatars.githubusercontent.com/u/78064891?v=4?s=100" width="100px;" alt="Yue Zhou"/><br /><sub><b>Yue Zhou</b></sub></a><br /><a href="https://github.com/NethermindEth/condor-rs/commits?author=Yue-Zhou1" title="Code">💻</a></td>
      <td align="center" valign="top" width="14.28%"><a href="https://github.com/NiDimi"><img src="https://avatars.githubusercontent.com/u/81875532?v=4?s=100" width="100px;" alt="Nick Dimitriou"/><br /><sub><b>Nick Dimitriou</b></sub></a><br /><a href="https://github.com/NethermindEth/condor-rs/pulls?q=is%3Apr+reviewed-by%3ANiDimi" title="Reviewed Pull Requests">👀</a></td>
    </tr>
  </tbody>
  <tfoot>
    <tr>
      <td align="center" size="13px" colspan="7">
        <img src="https://raw.githubusercontent.com/all-contributors/all-contributors-cli/1b8533af435da9854653492b1327a23a4dbd0a10/assets/logo-small.svg">
          <a href="https://all-contributors.js.org/docs/en/bot/usage">Add your contributions</a>
        </img>
      </td>
    </tr>
  </tfoot>
</table>

<!-- markdownlint-restore -->
<!-- prettier-ignore-end -->

<!-- ALL-CONTRIBUTORS-LIST:END -->

This project follows the [all-contributors](https://github.com/all-contributors/all-contributors) specification. Contributions of any kind welcome!