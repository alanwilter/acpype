# encoding: utf-8

__doc__ = """

Utility functions.

- directory handling
- pathname mangling
- running shell commands
- timer
- file checking
"""

import os
import stat
import tempfile
import copy
import re
import glob
import json
import subprocess
import time
import math
import pprint
import sys


def re_glob(dir_tag, reg_exp=""):
  fnames = glob.glob(dir_tag)
  return [f for f in fnames if re.search(reg_exp, f)]


def check_dir(dirname):
  if dirname.strip() == '':
    return
  if not os.path.isdir(dirname):
    os.makedirs(dirname)
    
  
def goto_dir(new_dir):
  if new_dir.strip() == '':
    return
  if not os.path.isdir(new_dir):
    os.makedirs(new_dir)
  os.chdir(new_dir)


def relpath(path):
  if path == '':
    return ''
  dirpath, basename = os.path.split(path)
  dirpath = os.path.relpath(dirpath, os.getcwd())
  if dirpath == '' or dirpath == '.':
    return basename
  return os.path.join(dirpath, basename)


def temp_fname(suffix=''):
  fd, fname = tempfile.mkstemp(suffix, 'tmp-', '.')
  f = os.fdopen(fd, 'w')
  f.close()
  os.unlink(fname)
  return os.path.basename(fname)


def fname_variant(fname):
  root, ext = os.path.splitext(fname)
  i = 1
  new_fname = "%s-%d%s" % (root, i, ext)
  while os.path.isfile(new_fname):
    i += 1
    new_fname = "%s-%d%s" % (root, i, ext)
  return new_fname


def clean_fname(*fnames):
  for fname in fnames:
    if os.path.isdir(fname):
      for root, dirs, files in os.walk(fname, topdown=False):
        for name in files:
          os.remove(os.path.join(root, name))
        for name in dirs:
          os.rmdir(os.path.join(root, name)) 
      os.rmdir(fname)
    elif os.path.isfile(fname):
      try:
        os.remove(fname)
      except:
        pass


def get_floats_from_string(s):
  val_strs = re.finditer(r'[-+]?([0-9]*\.[0-9]+|[0-9]+)', s)
  return [float(v.group()) for v in val_strs]
  
  
def write_dict(fname, d, indent=2):
  pprint.pprint(d, indent=2, stream=open(fname, "w"))


def read_dict(fname):
  try:
    txt = open(fname).read()
    d = eval(txt)
    if not isinstance(d, dict):
      raise Exception('Not a dictionary in ' + fname)
    return d
  except:
    raise Exception('Couldn\'t parse dictionary in ' + fname)
  

def is_same_dict_in_file(d, fname):
  try:
    saved_d = read_dict(fname)
    return saved_d == d
  except:
    return False


def words_in_file(fname):
  result = []
  for line in open(fname).readlines():
    result.extend(line.split())
  return result


def elapsed_time_str(time):
  s = str(time) + ' '
  minute = math.floor(time / 60.0)
  if minute > 60:
    hour = math.floor(minute / 60.0)
    partial_minute = math.fmod(time, 60.0)
    s += "%.f:%02.f:" % (hour, partial_minute)
  elif minute >= 1:
    s += "%.f:" % minute
  sec = math.fmod(time, 60.0)
  if sec < 0.01:
    s += "%07.4fs" % sec
  else:
    s += "%05.2fs" % sec
  return s
  
  
class Timer:
  def __init__(self):
    self._elapsed = 0;
    self._start = time.time()

  def start(self):
    self._start = time.time()
    self._elapsed = 0

  def stop(self):
    self._elapsed = time.time() - self._start

  def elapsed(self):
    if self._elapsed == 0:
      return time.time() - self._start
    else:
      return self._elapsed

  def str(self):
    elapsed_time = self.elapsed()
    return elapsed_time_str(elapsed_time)
    
  def __str__(self):
    return self.str()


def val_range(start, end, step):
  vals = []
  v = start
  while v <= end:
    vals.append(v)
    v += step
  return vals
 
 
class FileException(Exception):
  pass


def check_files(*fnames):
  """
  Checks for existence of fnames. Raises error if not found.
  """
  for fname in fnames:
    if not os.path.isfile(fname):
      raise FileException("Can't find {}".format(fname))


def which(program, other_binary_dirs=[]):
  """
  Reproduces Unix 'which' and looks in other_binary_dirs
  """

  def is_binary(fpath):
    return os.path.isfile(fpath) and os.access(fpath, os.X_OK)

  fpath, fname = os.path.split(program)
  if fpath:
    if is_binary(program):
      return program
  else:
    binary_dirs = os.environ["PATH"].split(os.pathsep)
    binary_dirs.extend(other_binary_dirs)
    for path in binary_dirs:
      exe_file = os.path.join(path, program)
      if is_binary(exe_file):
        return exe_file
  return None


def check_program(program):
  if not which(program):
    raise FileException("Can't find executable: " + program)


def check_output(fname, bad_words=[]):
  if not os.path.isfile(fname):
    raise FileException("Can't find output file: " + fname)
  txt = open(fname).read()
  if txt.strip() == '':
    raise Exception("Empty file: " + fname)
  for bad_word in bad_words:
    for i_line, line in enumerate(txt.splitlines()):
      if bad_word in line:
        raise Exception(
            "Output indicates %s error in line %d: %s" % (bad_word, i_line+1, fname))

def run_with_output(cmd):
  p = subprocess.Popen(
      cmd, 
      shell=True, 
      stdout=subprocess.PIPE, 
      stderr=subprocess.PIPE)
  return p.stdout.read()


def run_with_output_file(cmd, out_fname=None, in_fname=None):
  in_f = None
  out_f = None

  if in_fname and os.path.isfile(in_fname):
    in_f = open(in_fname)

  if out_fname:
    log_file = out_fname + '.log'
    out_f = open(log_file, 'w')
    sh_file = out_fname + '.sh'
    sh_cmd = cmd
    if in_f:
      sh_cmd += ' < ' + in_fname
    if log_file:
      sh_cmd += ' &> ' + log_file
    open(sh_file, 'w').write(sh_cmd)
    os.chmod(sh_file, stat.S_IRUSR | stat.S_IWUSR | stat.S_IXUSR)
    stopwatch = Timer()

  subprocess.call(
    cmd, 
    shell=True,
    stdin=in_f,
    stdout=out_f,
    stderr=out_f)

  if out_fname:
    stopwatch.stop()
    open(out_fname + '.time', 'w').write(stopwatch.str())
    out_f.close()



