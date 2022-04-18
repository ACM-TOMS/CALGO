# Contributing                                                                  {#page_contributing}

The RIDC software is managed by the GNU build system.  As such, the
developer release requires GNU autoconf, automake, libtool, m4, make
and their respective prerequisites.  If there are version mismatches
between the RIDC software and the local system, issuing the commands
`autoreconf -f` and `automake -a -c` should resolve version
errors and warning.  To build the documentation, Doxygen must be
installed, as well as appropriate Doxygen pre-requisites.  For
example, to build a PDF manual documenting the source code, Doxygen
requires a LaTeX compiler.

## Branching 

Contributors should fork
the git repository hosted at https://github.com/ongbw/ridc.git

If this project gets large enough, we will utilize
the [git-flow] workflow

Branch Name Pattern | Description
--------------------|------------
`master`            | tip of the `master` branch is always the latest stable release
`development`       | tip of the `development` branch is the current state of development and not expected to be stable or even usable
`feature/*`         | various feature branches are used to implement new features and should be based off the `development` branch
`release/*`         | a release branch is created from the `development` branch and used to prepare a new release and will be merged into `master`
`hotfix/*`          | hotfix branches are based off `master` or `development` to fix important and severe bugs and should be merged into `development` and `master` as soon as possible

Releases and release candidates are tagged in the form
`release-X.Y.Z(-RCa)`, where `X`, `Y`, and `Z` specify the version
with respect to [semantic versioning] and `a` the number of the
release candidate of that version.


## Commit Messages

Please keep commit messages clean and descriptive as possible.  The following
are suggested:

* Commit Title must not be longer than 50 characters

  If applicable, the title should start with a category name (such as
  `docu`, `tests`, ...)  followed by a colon (e.g. `"docu: add usage
  examples for RIDC"` ).

* Commit Description must have line wraps at 72 characters

* Please *sign* your commits (i.e. use `git commit -s`)



## How to Implement a New Feature?

-# create a fork/clone
-# switch to the `development` branch and pull in the latest changes
-# create a new branch `feature/XYZ` where `XYZ` is a short title of
   your planned feature (word seperation should be done with
   underscores, e.g. `feature/my_awesome_feature`)
-# hack and write Unit Tests
-# commit
-# repeat steps 4 and 5 until you feel your feature is in an almost
   usable state and most of the unit tests pass
-# write documentation for your feature
-# push your feature branch
-# stay tuned on reviews, remarks and suggestions by the other
   developers


[git-flow]: http://nvie.com/posts/a-successful-git-branching-model/

