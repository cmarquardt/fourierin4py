# -*- coding: utf-8 -*-

"""Calculate the current package version number based on git tags.

If possible, use the output of “git describe” modified to conform to the
versioning scheme that setuptools uses (see PEP 386).  Releases must be
labelled with annotated tags (signed tags are annotated) of the following
format:

    [v]<num>(.<num>)*[{a|b|c|rc}<num>(.<num>)*]

If “git describe” returns an error (likely because we're in an unpacked copy
of a release tarball, rather than a git repository), or returns a tag that
does not match the above format, the version is read from a RELEASE-VERSION
file.

To use this script, simply import it in your setup.py file and use the results
of get_version() as your package version:

    import version
    setup(
        version=version.get_version(),
        .
        .
        .
    )

This will automatically update the RELEASE-VERSION file.  The RELEASE-VERSION
file should *not* be checked into git but it *should* be included in sdist
tarballs (as should the version.py file).  To do this, run:

    echo include RELEASE-VERSION version.py >> MANIFEST.in
    echo RELEASE-VERSION >> .gitignore

With that setup, a new release can be labelled by simply invoking:

    git tag -s v1.0

The original version of this module is due to Douglas Creager. A nice 
write-up decribing the original script can be found at:

    http://dcreager.net/2010/02/10/setuptools-git-version-numbers/
    
This version includes several improvements made by Michal Nazarewicz (adaption
to PEP 386) and Ludwig Schwardt. I modified the output version string to be of 
the form

    <last tag>.post<N>+g<sha1>
    
which is exploiting the localised extension of the version strings as described
in PEP 400. <N> denotes the number of commits since the last tag. If no previous
tag is found, the script now returns '0.0' instead of aborting with an error. I 
also added support for tags named like '1.1' instead of 'v1.1', and changed the 
'git describe' options to allow for simple (not annotated) tags, as they seem to
occur in some of my repositories bridged to a svn repository using subgit.

"""

__author__ = ('Douglas Creager <dcreager@dcreager.net>',
              'Michal Nazarewicz <mina86@mina86.com>',
              'Ludwig Schwardt <ludwig@ska.ac.za>',
              'Christian Marquardt <christian@marquardt.sc>')
__license__ = 'This file is placed into the public domain.'
__maintainer__ = 'Christian Marquardt'
__email__ = 'christian@marquardt.sc'

__all__ = ('get_version')


import re
import subprocess
import sys


RELEASE_VERSION_FILE = 'RELEASE-VERSION'

# http://www.python.org/dev/peps/pep-0386/
_PEP386_SHORT_VERSION_RE = r'\d+(?:\.\d+)+(?:(?:[abc]|rc)\d+(?:\.\d+)*)?'
_PEP386_VERSION_RE = r'^%s(?:\.post\d+)?(?:\.dev\d+)?$' % (
    _PEP386_SHORT_VERSION_RE)
_GIT_DESCRIPTION_RE = r'^v{0,1}(?P<ver>%s)-(?P<commits>\d+)-g(?P<sha>[\da-f]+)$' % (
    _PEP386_SHORT_VERSION_RE)


def read_git_version():
    try:
        proc = subprocess.Popen(('git', 'describe', '--long', '--tags'),
                                stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        data, _ = proc.communicate()
        if proc.returncode:
            return
        ver = data.splitlines()[0].strip()
    except:
        return

    if not ver:
        return
    m = re.search(_GIT_DESCRIPTION_RE, ver)
    if not m:
        sys.stderr.write('version: git description (%s) is invalid, '
                         'ignoring\n' % (ver,))
        return

    commits = int(m.group('commits'))
    if not commits:
        return m.group('ver')
    else:
        return '%s.post%d+g%s' % (
            m.group('ver'), commits, m.group('sha'))


def read_release_version():
    try:
        fd = open(RELEASE_VERSION_FILE)
        try:
            ver = fd.readline().strip()
        finally:
            fd.close()
        if not re.search(_PEP386_VERSION_RE, ver):
            sys.stderr.write('version: release version (%s) is invalid, '
                             'will use it anyway\n' % (ver,))
        return ver
    except:
        return


def write_release_version(version):
    fd = open(RELEASE_VERSION_FILE, 'w+')
    fd.write('%s\n' % (version,))
    fd.close()


def get_version():
    release_version = read_release_version()
    version = read_git_version() or release_version
    if not version:
        version = '0.0'
    if version != release_version:
        write_release_version(version)
    return version


if __name__ == '__main__':
    print get_version()
