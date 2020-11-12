import tempfile
import shutil
import os
import itertools

class TemporaryDirectory:
    """Context Manager for working in a temporary directory"""

    def __init__(self, *args, **kwargs):
        self.temp_dir = tempfile.mkdtemp(*args, **kwargs)

    def __enter__(self):
        self.orig_dir = os.getcwd()
        os.chdir(self.temp_dir)
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        os.chdir(self.orig_dir)
        # If there was an error, do not delete the temporary
        # directory, so that the user can examine its contents
        if exc_type is None:
            shutil.rmtree(self.temp_dir)

    def copy_in(self, src, dst=None):
        """Copy a file into the temporary directory
        Arguments:
        src -- Source file name (relative to the original working directory)
        dst -- Destination file name (relative to the temporary directory)
               If not present, same as the source file name
        """
        if dst is None:
            dst = os.path.basename(src)
        if os.path.isabs(dst):
            raise ValueError("Destination path should not be absolute")
        abs_src = os.path.join(self.orig_dir, src)
        abs_dst = os.path.join(self.temp_dir, dst)
        shutil.copy(abs_src, abs_dst)
        return abs_dst

    def move_out(self, src, dst=None):
        """Move a file out of the temporary directory
        Arguments:
        src -- Source file name (relative to the temporary directory)
        dst -- Destination file name (relative to the original directory)
               If not present, same as the source file name
        """
        if os.path.isabs(src):
            raise ValueError("Source path should not be absolute")
        if dst is None:
            dst = src
        abs_src = os.path.join(self.temp_dir, src)
        abs_dst = os.path.join(self.orig_dir, dst)
        shutil.move(abs_src, abs_dst)
        return abs_dst

    def move_out_numbered(self, src, prefix, suffix):
        """Move a file out of the temporary directory, without overwriting old files
        Chooses a new file name based on the given prefix and suffix and a unique number
        Arguments:
        src -- Source file name (relative to the temporary directory)
        prefix -- prefix for the destination filename
        suffix -- suffix for the destination filename
        """
        if os.path.isabs(src):
            raise ValueError("Source path should not be absolute")
        if suffix.startswith('.'):
            suffix = suffix[1:]
        abs_src = os.path.join(self.temp_dir, src)
        abs_prefix = os.path.join(self.orig_dir, prefix)
        abs_dst = "%s.%s" % (abs_prefix, suffix)
        if os.path.exists(abs_dst):
            for i in itertools.count(start=1):
                abs_dst = "%s_%d.%s" % (abs_prefix, i, suffix)
                if not os.path.exists(abs_dst):
                    break
        shutil.move(abs_src, abs_dst)
        return abs_dst

