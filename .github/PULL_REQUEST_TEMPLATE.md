Thank you for contributing to `MoMA`!

Before submitting a pull request, please check the following:

- [ ] Does your pull request have tests? All pull requests are expected to
      increase (or at least not decrease) test coverage.
- [ ] Does your pull request pass `R CMD check` (or `devtools::check`) locally?
- [ ] If you're contributing a new feature, is it documented?
- [ ] Does every commit build and pass tests?
- [ ] Do your commit messages adhere to standard formatting?
      (Modified from the [Subsurface README](https://github.com/Subsurface-divelog/subsurface/blob/master/README.md).)

      Header line: explain the commit in one line (use the imperative)

      - Body of commit message is a few lines of text, explaining things
        in more detail, possibly giving some background about the issue
        being fixed, etc etc. Each "thought" should be in a separate
        bullet point (use hyphens for bullets and indent the body two
        spaces, like this text).

      - The body of the commit message can be several paragraphs, and
        please do proper word-wrap and keep columns shorter than about
        74 characters or so. That way "git log" will show things
        nicely even when it's indented.

      - Make sure you explain your solution and why you're doing what you're
        doing, as opposed to describing what you're doing. Reviewers and your
        future self can read the patch, but might not understand why a
        particular solution was implemented.

      - If your pull request is in response to an open issue, don't close it
        in your commit messages. We will close it in the merge commit message.

After submitting your pull request, please acknowledge
(in an comment on the pull request) that you have the right
to submit your code and you agree to license it under the GPL v2 or later.

If this is your first time contributing to the `MoMA` package, please also
add yourself to the
[`CONTRIBUTORS.md`](https://github.com/michaelweylandt/MoMA/blob/master/CONTRIBUTORS.md)
file.
