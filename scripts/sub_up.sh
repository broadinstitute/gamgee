# script to ensure that submodules reflect the latest versions of their respective .gitmodules branches
# affects the working directory only.

# warning: will remove any local changes to submodules

git submodule foreach --recursive git fetch --prune
git submodule sync
git submodule update --remote --recursive --force
