---
layout: post
title: "htslib tracking policy"
description: "how to interact with htslib"
category: setup
tags: [setup, submodule, best-practice]
---
{% include JB/setup %}

For now we are tracking their develop branch. The reason for this was because they hadn't updated their master branch in over a year and basically it had no BCF or CRAM support. Now that they have merged and promised to make more frequent releases, we could reconsider changing to master.

I personally think that following develop is a better option because it will give us access to the latest features. As long as we keep really good tests in gamgee, we will catch inconsistencies that they won't and we will be able to contribute them back. Their test infrastructure is pretty miniscule. If we choose to stay on master, we will risk falling behind significantly, not participate in the evolution of the library and the formats, and have to make huge modifications/adaptations every time we merge a big push from develop->master. So my vote is to **keep tracking develop**.

Regardless of what we track, we need certain discipline to work in this fashion. Here I lay out the procedures on how to work on our htslib repo.

Branch organization
---------------------------

### Master branch
Pretty much unused, since we can't force-push to it. This branch will be a simple mirror of samtools/htslib/master.

### Broad branch
This is our branch that we include as a submodule in gamgee. This will include their **develop** or **master** (whatever we decide), plus our modifications that *may or may not* be awaiting to be incorporated into samtools/htslib. This branch **will be force-pushed to** to keep our changes on top of their changes every time we rebase it against master or develop (whatever is the chosen base)

Making Changes
----------------------------
### To our repository (broad branch)
- Create a branch from broad
- make a pull-request against origin/broad

### Proposing changes to samtools/htslib
After the changes are accepted in our repo, you have to make a new branch to make a pull-request at htslib. We can't just pull-request the broad branch because other people's pull-requests to broad will show up in our pull-request to samtools/htslib. Therefore the process involves:
- create a new branch from samtools/htslib/develop with the name of your feature
- cherry pick your commits (from origin/broad) on top of it
- push to origin
- make a pull-request against samtools/htslib/develop
- (*if we are tracking samtools/htslib/develop*) after your pull-request is accepted, it is your responsibility to rebase origin/broad against the new version of their develop branch, so our modifications are always on top. Resolve any conflicts.
- update gamgee accordingly.
