import multiprocessing
from multiprocessing import Process, Pipe
from itertools import izip


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
