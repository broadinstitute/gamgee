echo -e "Publishing doxygen...\n"
git config --global user.email "travis@travis-ci.org"
git config --global user.name "travis-ci"
git clone --branch=gh-pages https://${GH_TOKEN}@github.com/broadinstitute/gamgee gh-pages 
cd gh-pages
mv $HOME/dox/html doxygen/
git commit -am "Latest doxygen documentation on successful travis build $TRAVIS_BUILD_NUMBER auto-pushed"
git push origin gh-pages
echo -e "Published doxygen.\n"
