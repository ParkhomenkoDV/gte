import matplotlib.pyplot as plt


def getValue(x1, x2, t):
    return x1 + (x2 - x1) * t


def Input(x, y):
    pointNum = int(input("Please enter the number of control points:"))
    for i in range(pointNum):
        x.append(float(input("Please enter the first" + str(i + 1) + "The x value of a control point:")))
        y.append(float(input("Please enter the first" + str(i + 1) + "Y value of a control point:")))
    return


def PointTest():
    for i in range(len(x)):
        plt.text(x[i], y[i], 'P' + str(i))
    return


def BezierPoint():
    for i in range(len(xLast)):
        plt.plot(xLast[i], yLast[i], '*')
    return


def Draw(x, y, num):
    if num == 3:
        x0 = getValue(x[0], x[1], i / cnt)
        y0 = getValue(y[0], y[1], i / cnt)
        x1 = getValue(x[1], x[2], i / cnt)
        y1 = getValue(y[1], y[2], i / cnt)
        plt.plot([x0, x1], [y0, y1])
        xLast.append(getValue(x0, x1, i / cnt))
        yLast.append(getValue(y0, y1, i / cnt))
        plt.plot(xLast[i], yLast[i], '*')
        plt.pause(0.1)
    else:
        xNext = []
        yNext = []
        for j in range(num - 1):
            xNext.append(getValue(x[j], x[j + 1], i / cnt))
            yNext.append(getValue(y[j], y[j + 1], i / cnt))
        plt.plot(xNext, yNext)
        plt.pause(0.1)
        Draw(xNext, yNext, num - 1)
    return


def Show():
    plt.grid(True)
    plt.plot(x, y)
    PointTest()
    BezierPoint()
    return


#
x = []
y = []
xLast = []
yLast = []
cnt = 10
Input(x, y)
fig = plt.figure()
ax = fig.add_subplot(111)
Show()
for i in range(cnt + 1):
    Draw(x, y, len(x))
    plt.pause(0.1)
    plt.cla()
    Show()
plt.cla()
Show()
plt.plot(xLast, yLast)
plt.show()
