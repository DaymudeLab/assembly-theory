import multiprocessing


def wrapper(queue, func, args, kwargs):
    """Executes the function and puts the result in the queue."""
    try:
        result = func(*args, **kwargs)
        queue.put(result)
    except Exception as e:
        queue.put(e)


def run_with_timeout(func, timeout, *args, **kwargs):
    """Runs a function with a timeout, raising an exception if it exceeds the limit."""
    queue = multiprocessing.Queue()
    process = multiprocessing.Process(target=wrapper, args=(queue, func, args, kwargs))
    process.start()
    process.join(timeout)

    if process.is_alive():
        process.terminate()
        process.join()
        raise TimeoutError(
            f"Function {func.__name__} exceeded time limit of {timeout} seconds."
        )

    if not queue.empty():
        result = queue.get()
        if isinstance(result, Exception):
            raise result  # Re-raise any function exceptions
        return result

    raise TimeoutError("Function did not return any result.")
