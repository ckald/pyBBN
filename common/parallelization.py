import multiprocessing
from multiprocessing import Process, Pipe, Pool
from itertools import izip

target_func = None


def initializer(target_func):
    global target
    target = target_func


def target(args, **kw):
    object, method_name = args[0], args[1]
    args = args[2:]
    result = getattr(object, method_name)(*args)
    return result

pool = Pool(processes=7, initializer=initializer, initargs=(target,))


def spawn(f):
    def fun(pipe, x):
        pipe.send(f(x))
        pipe.close()
    return fun


def parmap(f, X, workers=multiprocessing.cpu_count()-1):
    pipe = [Pipe() for x in X]
    processes = [Process(target=spawn(f), args=(c, x)) for x, (p, c) in izip(X, pipe)]
    numProcesses = len(processes)
    processNum = 0
    outputList = []
    while processNum < numProcesses:
        endProcessNum = min(processNum+workers, numProcesses)
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

    result = pool.map(func, arguments)
    return result
