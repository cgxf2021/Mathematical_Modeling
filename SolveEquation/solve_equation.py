class SolveEquation():
    """ A class for solving equation """
    
    def __init__(self):
        self.xc = 0.0

    
    def bisect(self, fun, a, b, tol):
        """ 
            A method of solving the roots of an equation by dichotomy

            fun: Target equation
            a: Left end of interval
            b: Right end of interval
            c: The midpoint of the interval
            tol: tolerance (0.00005)
            return self.xc: Approximate solution
        """

        try:
            if fun(a) * fun(b) >= 0:
                raise Exception()
        except Exception:
            print('f(a)*f(b)<0 not saitisfied')
            return None
        while (b - a) / 2 > tol:
            c = (a + b) / 2
            if fun(c) == 0:
                break
            else:
                if fun(a) * fun(c) < 0:
                    b = c
                else:
                    a = c
        self.xc = (a + b) / 2
        return self.xc


    def fpi(self, fun, x0, k):
        """
            A method of solving the roots of an equation by iteration of fixed points

            fun: Target equation (f(x) = x)
            x0: Initial estimate
            k: The number of iterations
            return self.xc: Approximate solution
        """
        x = x0
        for i in range(k):
            x = fun(x)
        self.xc = x
        return self.xc


    def newton(self, fun, dirf, x0, k, m=1):
        """
            A method of solving the roots of an equation by Newton's method

            fun: Target equation [f(x)]
            dirf: The derivative of target equation [f'(x)]
            x0: Initial estimate
            k: The number of iterations
            m: m multiple roots
            return self.xc: Approximate solution
        """
        x = x0
        for i in range(k):
            try:
                x = x - m * fun(x) / dirf(x)
            except ZeroDivisionError:
                print('The denominator is not 0')
                return x
            else:
                if fun(x) == 0:
                    break
        self.xc = x
        return self.xc
    
    
    def secant(self, fun, x0, x1, k):
        """
            A method of solving the roots of an equation by secant method

            fun: Target equation [f(x)]
            x0: One initial estimate
            x1: The other initial estimate
            k: The number of iterations
            return self.xc: Approximate solution
        """
        for i in range(k):
            try:
                x = x1 - fun(x1) * (x1 - x0) / (fun(x1) - fun(x0))
            except ZeroDivisionError:
                print('The denominator is not 0')
                return x
            else:
                if fun(x) == 0:
                    break
                x0 = x1
                x1 = x
        self.xc = x
        return self.xc


    def regula(self, fun, a, b, k):
        """
            A method of solving the roots of an equation by secant method

            fun: Target equation [f(x)]
            a: The first initial estimate
            b: The second initial estimate
            k: The number of iterations
            retuen self.xc: Approximate solution
        """
        try:
            if fun(a) * fun(b) >= 0:
                raise Exception()
        except Exception:
            print('f(a)*f(b)<0 not saitisfied')
            return None
        for i in range(k):
            try:
                c = (b * fun(a) - a * fun(b)) / (fun(a) - fun(b))
            except ZeroDivisionError:
                print('The denominator is not 0')
                return c
            else:
                if fun(c) == 0:
                    break
                if fun(a) * fun(c) < 0:
                    b = c
                else:
                    a = c
        self.xc = c
        return self.xc
            