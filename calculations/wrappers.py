import colorama

colorama.init(autoreset=True)

def running(function):
    def wrapper():
        print(colorama.Fore.YELLOW + f'Running {func.__name__}')
        func()
        # print('Complete')
    return wrapper

def tictoc(function):
    def wrapper():
        t1 = time.time()
        function()
        t2 = time.time() - t1
        print('{0} ran in {1} seconds'.format(function.__name__, t2))

    return wrapper