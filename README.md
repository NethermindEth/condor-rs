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
      <td align="center" valign="top" width="14.28%"><a href="https://github.com/JamesEBall"><img src="https://avatars.githubusercontent.com/u/73405337?v=4?s=100" width="100px;" alt="James Ball"/><br /><sub><b>James Ball</b></sub></a><br /><a href="#ideas-JamesEBall" title="Ideas, Planning, & Feedback">🤔</a> <a href="#fundingFinding-JamesEBall" title="Funding Finding">🔍</a> <a href="#research-JamesEBall" title="Research">🔬</a></td>
      <td align="center" valign="top" width="14.28%"><a href="https://github.com/omibo"><img src="https://avatars.githubusercontent.com/u/42227752?v=4?s=100" width="100px;" alt="Omid Bodaghi"/><br /><sub><b>Omid Bodaghi</b></sub></a><br /><a href="https://github.com/NethermindEth/condor-rs/pulls?q=is%3Apr+reviewed-by%3Aomibo" title="Reviewed Pull Requests">👀</a></td>
    </tr>
    <tr>
      <td align="center" valign="top" width="14.28%"><a href="https://github.com/v0-e"><img src="https://avatars.githubusercontent.com/u/134806759?v=4?s=100" width="100px;" alt="v0-e"/><br /><sub><b>v0-e</b></sub></a><br /><a href="https://github.com/NethermindEth/condor-rs/pulls?q=is%3Apr+reviewed-by%3Av0-e" title="Reviewed Pull Requests">👀</a></td>
      <td align="center" valign="top" width="14.28%"><a href="https://maksimryndin.github.io/"><img src="https://avatars.githubusercontent.com/u/16288656?v=4?s=100" width="100px;" alt="maksimryndin"/><br /><sub><b>maksimryndin</b></sub></a><br /><a href="https://github.com/NethermindEth/condor-rs/pulls?q=is%3Apr+reviewed-by%3Amaksimryndin" title="Reviewed Pull Requests">👀</a></td>
      <td align="center" valign="top" width="14.28%"><a href="https://github.com/earthling1984"><img src="https://avatars.githubusercontent.com/u/19665196?v=4?s=100" width="100px;" alt="Mithun Vaidhyanathan"/><br /><sub><b>Mithun Vaidhyanathan</b></sub></a><br /><a href="#research-earthling1984" title="Research">🔬</a></td>
      <td align="center" valign="top" width="14.28%"><a href="https://github.com/smartprogrammer93"><img src="https://avatars.githubusercontent.com/u/33181301?v=4?s=100" width="100px;" alt="Ahmad Bitar"/><br /><sub><b>Ahmad Bitar</b></sub></a><br /><a href="#ideas-smartprogrammer93" title="Ideas, Planning, & Feedback">🤔</a></td>
      <td align="center" valign="top" width="14.28%"><a href="https://github.com/HendrikWaldner"><img src="https://avatars.githubusercontent.com/u/33893964?v=4?s=100" width="100px;" alt="HendrikWaldner"/><br /><sub><b>HendrikWaldner</b></sub></a><br /><a href="#research-HendrikWaldner" title="Research">🔬</a></td>
      <td align="center" valign="top" width="14.28%"><a href="https://www.linkedin.com/in/aragirtas"><img src="https://avatars.githubusercontent.com/u/95563727?v=4?s=100" width="100px;" alt="Ahmet Ramazan Agirtas"/><br /><sub><b>Ahmet Ramazan Agirtas</b></sub></a><br /><a href="#research-ahmetramazan" title="Research">🔬</a></td>
      <td align="center" valign="top" width="14.28%"><a href="https://github.com/mpzajac"><img src="https://avatars.githubusercontent.com/u/25900405?v=4?s=100" width="100px;" alt="mpzajac"/><br /><sub><b>mpzajac</b></sub></a><br /><a href="#research-mpzajac" title="Research">🔬</a></td>
    </tr>
    <tr>
      <td align="center" valign="top" width="14.28%"><a href="https://matthewklein.co.uk/"><img src="https://avatars.githubusercontent.com/u/96837318?v=4?s=100" width="100px;" alt="Matthew Klein"/><br /><sub><b>Matthew Klein</b></sub></a><br /><a href="#research-matthew-a-klein" title="Research">🔬</a></td>
      <td align="center" valign="top" width="14.28%"><a href="https://github.com/GustaveCharles"><img src="https://avatars.githubusercontent.com/u/91184289?v=4?s=100" width="100px;" alt="GustaveCharles"/><br /><sub><b>GustaveCharles</b></sub></a><br /><a href="#research-GustaveCharles" title="Research">🔬</a></td>
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