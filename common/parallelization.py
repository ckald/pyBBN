from itertools import izip
from multiprocessing import Pool, Pipe, Process


def target(args, **kw):
    obj, method_name = args[0], args[1]
    args = args[2:]
    result = getattr(obj, method_name)(*args, **kw)
    return result

worker_count = 50  # multiprocessing.cpu_count() / 2
pool = None
orders = []


def init_pool():
    global pool
    pool = Pool(processes=worker_count)


def close_pool():
    pool.close()
    pool.terminate()
    pool.join()


def map_orders():
    results = []
    for i, (particle, method, arg) in enumerate(orders):
        results.append((i, parmap(lambda x: getattr(particle, method)(x), arg)))
    return results


def poolmap(cls, func_name, arguments):
    arguments = [(cls, func_name, arg) for arg in arguments]
    func = target

    result = pool.map_async(func, arguments)
    return result


def spawn(f):
    def fun(pipe, x):
        pipe.send(f(x))
        pipe.close()
    return fun


def parmap(f, X, workers=worker_count):
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
