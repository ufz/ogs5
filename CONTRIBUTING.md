# How to contribute



## Getting Started

* Make sure you have a [GitHub account](https://github.com/signup/free)
* Fork the repository on GitHub
* Install clang-format (download from [here](http://llvm.org/builds/) for Windows)

## Making Changes

* Create a topic branch from where you want to base your work.
  * This is usually the master branch.
  * Only target other branches if you are certain your fix must be on that
    branch.
  * To quickly create a topic branch based on master; `git branch
    my_topic_branch master` then checkout the new branch with `git
    checkout my_topic_branch`.  Please avoid working directly on the
    `master` branch.
* Make commits of logical units.
* Make sure your code conforms to the [styleguide][styleguide].
* Format the codes using `clang-format -i -style=file <file name>` before committing. To format all the files, you can use scripts/clang-format/run-clang-format.sh (or .bat for Windows).


## Submitting Changes

* Push your changes to a topic branch in your fork of the repository.
* Submit a pull request to the main repository.

# Additional Resources

* [General GitHub documentation](http://help.github.com/)
* [GitHub pull request documentation](http://help.github.com/send-pull-requests/)
* [OGS Styleguide][styleguide]

[styleguide]: http://ufz.github.com/styleguide/cppguide.xml