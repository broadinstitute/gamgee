---
layout: page
title: Gamgee
tagline: native sequencing data processing
---
{% include JB/setup %}

[![Build Status](https://travis-ci.org/broadinstitute/gamgee.svg?branch=master)](https://travis-ci.org/broadinstitute/gamgee)
[![Coverage Status](https://coveralls.io/repos/broadinstitute/gamgee/badge.png)](https://coveralls.io/r/broadinstitute/gamgee)
[![Issue Stats](http://issuestats.com/github/broadinstitute/gamgee/badge/pr)](http://issuestats.com/github/broadinstitute/gamgee)
[![Issue Stats](http://issuestats.com/github/broadinstitute/gamgee/badge/issue)](http://issuestats.com/github/broadinstitute/gamgee)

### What is Gamgee

A C++ interface for key next generation sequencing file formats. This library
provides an API for handilng FASTA, FASTQ, SAM, BAM, CRAM, VCF, BCF, BED,
GATK/Picard intervals and interval_list files.

The APIs are free to change until we reach version 1.x

<a class="btn danger" href="http://broadinstitute.github.io/gamgee/doxygen/">Developer docs</a>
<a class="btn primary" href="http://github.com/MauricioCarneiro/gamgee/">View on Github</a>
<a class="btn info" href="https://github.com/MauricioCarneiro/gamgee/zipball/master">Download ZIP</a>
<a class="btn info" href="https://github.com/MauricioCarneiro/gamgee/tarball/master">Download TAR</a>

### Developer Blog

<ul class="posts">
  {% for post in site.posts %}
    <li><span>{{ post.date | date_to_string }}</span> &raquo; <a href="{{ BASE_PATH }}{{ post.url }}">{{ post.title }}</a></li>
  {% endfor %}
</ul>

### Using Gamgee

Instructions on how to set up your development environment, building, running and testing are [here](setup/20140924-setting-up-your-development-environment/)

### Licensing

Gamgee is copyrighted under the [MIT License](http://opensource.org/licenses/MIT).
