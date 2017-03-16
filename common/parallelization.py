from itertools import izip
from multiprocessing import Pool, Pipe, Process, cpu_count


worker_count = cpu_count() - 1
orders = []


# ## Pool-oriented parallelization

pool = None


def init_pool(workers=None):
    global pool
    global worker_count

    if workers is not None:
        worker_count = min(workers, worker_count)
    pool = Pool(processes=worker_count)


def close_pool():
    pool.close()
    pool.terminate()
    pool.join()


def target(args, **kw):
    obj, method_name = args[0], args[1]
    args = args[2:]
    result = getattr(obj, method_name)(*args, **kw)
    return result


def poolmap(cls, func_name, arguments):
    arguments = [(cls, func_name, arg) for arg in arguments]
    func = target

    result = pool.map_async(func, arguments)
    return result


# ## Worker-oriented parallelization

def map_order(i):
    particle, method, arg = orders[i]
    return getattr(particle, method)(arg)


def map_orders():
    return parmap(map_order, range(len(orders)))


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
