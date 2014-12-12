import sys
import time
import codecs
from functools import wraps
from collections import deque


class PicklableObject(object):

    _saveable_fields = None

    def __getstate__(self):
        if self._saveable_fields:
            return {key: value for key, value in self.__dict__.items()
                    if key in self._saveable_fields}
        return self.__dict__

    def __setstate__(self, data):
        if self._saveable_fields:
            data = {key: value for key, value in data.items()
                    if key in self._saveable_fields}
        self.__dict__.update(data)


class Logger(object):
    """ Convenient double logger that redirects `stdout` and save the output also to the file """
    def __init__(self, filename):
        self.terminal = sys.stdout
        self.log = codecs.open(filename, "w", encoding="utf8")

    def write(self, message, terminal=True, log=True):
        if terminal:
            self.terminal.write(message)
        if log:
            self.log.write(message.decode('utf8'))

    def __del__(self):
        if sys:
            sys.stdout = self.terminal

    def __getattr__(self, attr):
        return getattr(self.terminal, attr)


class ring_deque(deque):
    """ Circular deque implementation """

    max = 0
    min = 0

    def __init__(self, data, length):
        self.length = length
        super(ring_deque, self).__init__(data)

    def append_more(self, data):

        if len(self) > self.length:
            self.popleft()

        self.max = max(data, self.max)
        self.min = min(data, self.min)
        super(ring_deque, self).append(data)

    def append(self, data):
        self.max = data
        self.min = data

        self.append = self.append_more

        super(ring_deque, self).append(data)


class benchmark(object):
    """ Simple benchmarking context manager """
    def __init__(self, name):
        self.name = name

    def __enter__(self):
        self.start = time.time()

    def __exit__(self, ty, val, tb):
        end = time.time()
        print("%s : %0.5f seconds" % (self.name, end-self.start))
        return False


def echo(func):
    """ Print the return value of the function """
    @wraps(func)
    def wrapper(*args, **kw):
        val = func(*args, **kw)
        print val
        return val
    return wrapper


def memodict(f):
    """ Memoization decorator for a function taking a single argument """
    class memodict(dict):
        def __missing__(self, key):
            ret = self[key] = f(key)
            return ret
    return memodict().__getitem__


class MemoizeMutable:
    """ Multiple arguments implementation of the memoization decorator """
    def __init__(self, fn):
        self.fn = fn
        self.memo = {}

    def __call__(self, *args):
        pickle = tuple(args)
        if not pickle in self.memo:
            self.memo[pickle] = self.fn(*args)

        return self.memo[pickle]
