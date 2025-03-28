Dev team's guide to `assembly-theory`
***

The maintainers will govern the project using the committee model: high-level decisions about the project's direction require maintainer consensus, major code changes require majority approval, hotfixes and patches require just one maintainer approval, new maintainers can be added by unanimous decision of the existing maintainers, and existing maintainers can step down with advance notice.
After the initial development period, maintainers are expected to perform issue prioritization and code review on a reasonably frequent basis and meet at least quarterly to discuss larger project decisions.

All source code will remain openly available on GitHub throughout the lifetime of the project, supporting user feedback and discussions via GitHub Issues, change requests and code review via GitHub Pull Requests, and continuous integration and deployment via GitHub Actions.


- Main dev team: work in branches. Rebase your branch onto main to update your
  branch. Only merge from main if your code significantly deviates from main.
  Make a PR to merge changes to main.

- Your PR won't be accepted if it doesn't pass **all** clippy and rustfmt
  checks.

- You should read through the [Rust API guidelines
  checklist](https://rust-lang.github.io/api-guidelines/checklist.html).
  Your code should follow these guidelines closely.
