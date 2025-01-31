Dev team's guide to ORCA
***

- Main dev team: work in branches. Rebase your branch onto main to update your
  branch. Only merge from main if your code significantly deviates from main.
  Make a PR to merge changes to main.

- Your PR won't be accepted if it doesn't pass **all** clippy and rustfmt
  checks.

- You should read through the [Rust API guidelines
  checklist](https://rust-lang.github.io/api-guidelines/checklist.html).
  Your code should follow these guidelines closely.
