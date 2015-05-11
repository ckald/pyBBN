import multiprocessing
from multiprocessing import Process, Pipe, Pool
from itertools import izip


def target(args, **kw):
    obj, method_name = args[0], args[1]
    args = args[2:]
    result = getattr(obj, method_name)(*args, **kw)
    return result

worker_count = multiprocessing.cpu_count() / 2
pool = Pool(processes=worker_count)


def spawn(f):
    def fun(pipe, x):
        pipe.send(f(x))
        pipe.close()
    return fun


def parmap(f, X):
    pipe = [Pipe() for x in X]
    processes = [Process(target=spawn(f), args=(c, x)) for x, (_, c) in izip(X, pipe)]
    numProcesses = len(processes)
    processNum = 0
    outputList = []
    while processNum < numProcesses:
        endProcessNum = min(processNum+worker_count, numProcesses)
        for proc in processes[processNum:endProcessNum]:
            proc.start()
        for proc in processes[processNum:endProcessNum]:
            proc.join()
        for proc, c in pipe[processNum:endProcessNum]:
            outputList.append(proc.recv())
        processNum = endProcessNum
    return outputList


def poolmap(cls, func_name, arguments):
    arguments = [(cls, func_name, arg) for arg in arguments]
    func = target

    result = pool.map_async(func, arguments, chunksize=5)
    return result
