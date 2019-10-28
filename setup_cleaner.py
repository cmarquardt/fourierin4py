# setup_cleaner
# -------------
#
# Setup_cleaner is a localised version of janitor.py from https://github.com/maphew/setupext-janitor

from distutils import dir_util, errors
from distutils.command.clean import clean as _CleanCommand
import os.path
import traceback

debug = False


class CleanCommand(_CleanCommand):
    """
    Extend the clean command to do additional house keeping.

    The traditional distutils clean command removes the by-products of
    compiling extension code.  This class extends it to remove the
    similar by-products generated by producing a Python distribution.
    Most notably, it will remove .egg/.egg-info directories, the
    generated distribution, those pesky *__pycache__* directories,
    and even the virtual environment that it is running in.

    The level of cleanliness is controlled by command-line options as
    you might expect.  The targets that are removed are influenced by
    the commands that created them.  For example, if you set a custom
    distribution directory using the ``--dist-dir`` option or the
    matching snippet in *setup.cfg*, then this extension will honor
    that setting.  It even goes as far as to detect the virtual
    environment directory based on environment variables.

    This all sounds a little dangerous... there is little to worry
    about though.  This command only removes what it is configured to
    remove which is nothing by default.  It also honors the
    ``--dry-run`` global option so that there should be no question
    what it is going to remove.

    """

    # See _set_options for `user_options`

    def initialize_options(self):
        _CleanCommand.initialize_options(self)
        self.build = False
        self.dist = False
        self.eggs = False
        self.egg_base = None
        self.environment = False
        self.pycache = False
        self.pytest_cache = False
        self.virtualenv_dir = None

    def finalize_options(self):
        _CleanCommand.finalize_options(self)
        try:
            self.set_undefined_options(
                'egg_info', ('egg_base', 'egg_base'))
        except errors.DistutilsError:
            pass

        if self.egg_base is None:
            self.egg_base = os.curdir

        if self.all:
            for flag in self.boolean_options:
                setattr(self, flag, True)

        if self.environment and self.virtualenv_dir is None:
            self.virtualenv_dir = os.environ.get('VIRTUAL_ENV', None)

    def run(self):
        _CleanCommand.run(self)

        dir_names = set()

        if self.dist:
            for cmd_name, _ in self.distribution.get_command_list():
                if 'dist' in cmd_name:
                    command = self.distribution.get_command_obj(cmd_name)
                    # command.ensure_finalized()
                    # Stop premature exit for RPM-on-NT err
                    # https://github.com/dave-shawley/setupext-janitor/issues/12
                    try:
                        command.ensure_finalized()
                    except Exception as err:
                        skip = "don't know how to create RPM distributions on platform nt"
                        if skip in err.args:
                            print('-' * 50, '\nException encountered and ignored:')
                            print('{} {}'.format(err.__class__.__name__, err.args[0]))
                            if debug: traceback.print_exc()
                            print('-' * 50)
                        else:
                            raise err

                    if getattr(command, 'dist_dir', None):
                        dir_names.add(command.dist_dir)

        if self.build:
            for root, dirs, _ in os.walk(os.curdir):
                if 'build' in dirs:
                    dir_names.add(os.path.join(root, 'build'))

        if self.eggs:
            for name in os.listdir(self.egg_base):
                if name.endswith('.egg-info'):
                    dir_names.add(os.path.join(self.egg_base, name))
            for name in os.listdir(os.curdir):
                for e in ['.egg', '.eggs']:
                    if name.endswith(e):
                        dir_names.add(name)

        if self.environment and self.virtualenv_dir:
            dir_names.add(self.virtualenv_dir)

        if self.pycache:
            for root, dirs, _ in os.walk(os.curdir):
                if '__pycache__' in dirs:
                    dir_names.add(os.path.join(root, '__pycache__'))

        if self.pytest:
            for root, dirs, _ in os.walk(os.curdir):
                if '.pytest_cache' in dirs:
                    dir_names.add(os.path.join(root, '.pytest_cache'))

        for dir_name in dir_names:
            if os.path.exists(dir_name):
                dir_util.remove_tree(dir_name, dry_run = self.dry_run)
            else:
                self.announce(
                    'skipping {0} since it does not exist'.format(dir_name))

        # Remove .pyc files (always)

        for root, _, files in os.walk(os.curdir):
            for basename in files:
                for ext in ['.pyc', '.so']:
                    if basename.endswith(ext):
                        os.remove(os.path.join(root, basename))



def _set_options():
    """
    Set the options for CleanCommand.

    There are a number of reasons that this has to be done in an
    external function instead of inline in the class.  First of all,
    the setuptools machinery really wants the options to be defined
    in a class attribute - otherwise, the help command doesn't work
    so we need a class attribute.  However, we are extending an
    existing command and do not want to "monkey patch" over it so
    we need to define a *new* class attribute with the same name
    that contains a copy of the base class value.  This could be
    accomplished using some magic in ``__new__`` but I would much
    rather set the class attribute externally... it's just cleaner.

    """
    CleanCommand.user_options = _CleanCommand.user_options[:]
    CleanCommand.user_options.extend([
        ('build', 'B', 'remove build directory'),
        ('dist', 'd', 'remove distribution directory'),
        ('eggs', 'e', 'remove egg and egg-info directories'),
        ('pytest', None, 'remove .pytest_cache directories'),
        ('bytecode', None, 'remove python bytecode (__pycache__) directories'),
        ('environment', None, 'remove virtual environment directory'),

        ('egg-base=', None,
         'directory containing .egg-info directories '
         '(default: top of the source tree)'),
        ('virtualenv-dir=', None,
         'root directory for the virtual directory '
         '(default: value of VIRTUAL_ENV environment variable)'),
    ])
    CleanCommand.boolean_options = _CleanCommand.boolean_options[:]
    CleanCommand.boolean_options.extend(
        ['build', 'dist', 'eggs', 'pycache','pytest'])


_set_options()