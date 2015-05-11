import multiprocessing
from multiprocessing import Pool


def target(args, **kw):
    obj, method_name = args[0], args[1]
    args = args[2:]
    result = getattr(obj, method_name)(*args, **kw)
    return result

worker_count = multiprocessing.cpu_count() / 2
pool = Pool(processes=worker_count)


def poolmap(cls, func_name, arguments):
    arguments = [(cls, func_name, arg) for arg in arguments]
    func = target

    result = pool.map_async(func, arguments, chunksize=5)
    return result
