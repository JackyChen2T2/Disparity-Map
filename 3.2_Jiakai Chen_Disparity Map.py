# import MATLAB plotting library for visual I/O and processing
import matplotlib
matplotlib.use('Qt5Agg')
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.image as matimage
import matplotlib.pyplot as matplot
# import math library
import math
# for cloning matrix
from copy import copy, deepcopy


def ReadFile(i, window):
    """Function that reads txt and image files noted by a number i\n
    Return 2 matrices of cropped images"""
    # assign file names noted by i
    dataFileName = str(i) + "_calib.txt"
    leftImageName = str(i) + "_im0.png"
    rightImageName = str(i) + "_im1.png"
    
    # read calib.txt for photo sizes
    try:
        width = 0
        height = 0
        cam0 = 0
        cam1 = 0
        leftImage = []
        rightImage = []
        with open(dataFileName, "r") as file:
            for line in file:
                lineSeg = line.split("=")
                # extract displacement between 2 cameras
                if lineSeg[0] == "cam0":
                    temp = lineSeg[1].split(" ")
                    cam0 = math.floor(float(temp[2][:-1])/window) * window
                elif lineSeg[0] == "cam1":
                    temp = lineSeg[1].split(" ")
                    cam1 = math.floor(float(temp[2][:-1])/window) * window
                # round off last digits to ensure a common denominator of 10
                elif lineSeg[0] == "width":
                    width = math.floor(int(lineSeg[1])/window) * window
                elif lineSeg[0] == "height":
                    height = math.floor(int(lineSeg[1])/window) * window
        # read the cropped left and right images
        leftImage = (matimage.imread(leftImageName))[:height][:width]
        rightImage = (matimage.imread(rightImageName))[:height][:width]
        # calculate camera displacement as the largest possible displacement
        camDisplacement = cam1 - cam0
    except:
        print("Error in file reading number " + str(i))
    return [leftImage, rightImage, camDisplacement]


def Variance(image, window):
    """Function that calculates regional directional variances of an image\n
    Return an matrix of variances with a size scaled down by window size"""
    # assign image and matrix dimensions and preallocate the matrix
    imageWidth = len(image[0])
    imageHeight = len(image)
    matrixWidth = imageWidth // window
    matrixHeight = imageHeight // window
    matri = []
    # calculate regional variance of the image and store values to the matrix
    for m_row in range(matrixHeight):
        sub_matri = []
        for m_col in range(matrixWidth):
            i_left = m_col * window
            i_right = i_left + window - 1
            i_top = m_row * window
            i_bottom = i_top + window - 1
            rgb_matri = []
            worthy = False
            for c in range(3):
                I1 = 0
                I2 = 0
                I3 = 0
                I4 = 0
                for i_row in range(i_top+1, i_bottom, 1):
                    for i_col in range(i_left, i_right, 1):
                        I1 += (image[i_row][i_col][c] - image[i_row+1][i_col][c]) ** 2
                        I2 += (image[i_row][i_col][c] - image[i_row][i_col+1][c]) ** 2
                        I3 += (image[i_row][i_col][c] - image[i_row+1][i_col+1][c]) ** 2
                        I4 += (image[i_row][i_col][c] - image[i_row-1][i_col+1][c]) ** 2
                # modification: instead of detecting interesting points, detect edge points --> ease the color filling ***************
                temp = I1 / (window**2) + I2 / (window**2) + I3 / (window**2) + I4 / (window**2)
                # experimental threshold ***************
                if temp > 0.0002:      #0.0002
                    worthy = True
                rgb_matri += [temp]
            # if variances of all rgb channels are too low, take out the point
            if worthy:
                sub_matri += [rgb_matri]
            else:
                sub_matri += [[0,0,0]]
        matri += [sub_matri]
        print("Variance: row " + str(m_row) + "/" + str(matrixHeight - 1))
    print("COMPLETED.")
    return matri


def Correlation(leftImage, rightImage, leftVariance, rightVariance, window, maxBound):              # Experimental: no filling in this function ****
    """Function that calculates correlation coefficient between each pair of variance points\n
    the disparity between a pair of best match is stored in a disparity matrix\n
    Return an disparity matrix of correlated points in the size of variance matrix"""
    varianceWidth = len(leftVariance[0])
    varianceHeight = len(leftVariance)
    # Prepare a pair of two residual matrices for correlation coefficient calculations
    leftResidual = Residual(leftImage, leftVariance, window, varianceWidth, varianceHeight)
    rightResidual = Residual(rightImage, rightVariance, window, varianceWidth, varianceHeight)
    # Calculate and find best correlation coefficients in each row of two image
    correlation = []
    for c_row in range(varianceHeight):
        subcorrel = []
        for c_col in range(varianceWidth):
            if leftResidual[c_row][c_col] == []:            # Experimental: if this point is not considered, mark it with -1 for later ***********
                subcorrel += [-1]
                continue
            
            #correlationCoefficient = [] --> correlation matrix
            # This is the disparity of most correlated point found for each color channel
            BEST_MATCH = 0
            # This is the correlation coefficient of the above most correlated points
            MATCHED_CC = -1
            # This is the max search bound estimate from the displacement between cameras, assuming inversed image plane is between camera and object
            bound = c_col - maxBound // window
            if bound < -1: bound = -1
            # assume no point remains on the same location in both images, as such a point is infinitely far away
            for right2left in range(c_col-1, bound, -1):
                if rightResidual[c_row][right2left] == []:
                    #correlationCoefficient += [[]] --> correlation matrix
                    continue
                coefficient = [0,0,0]
                num = [0,0,0]
                den = [0,0,0]
                den_l = [0,0,0]
                den_r = [0,0,0]
                avg = 0
                for c in range(3):
                    remaining = 3
                    for n in range(window**2):
                        num[c] += leftResidual[c_row][c_col][n][c] * rightResidual[c_row][right2left][n][c]
                        den_l[c] += leftResidual[c_row][c_col][n][c] ** 2
                        den_r[c] += rightResidual[c_row][right2left][n][c] ** 2
                    den[c] = (den_l[c] * den_r[c]) ** (1/2)
                    if den[c] == 0:
                        remaining -= 1
                    else:
                        coefficient[c] = num[c] / den[c]
                    avg += coefficient[c]
                if remaining == 0:
                    continue
                avg = avg / remaining                                   # experimental: averaging coreelations from 3 channels ******
                if avg > MATCHED_CC:
                    MATCHED_CC = avg
                    BEST_MATCH = c_col - right2left
                #correlationCoefficient += [coefficient] --> correlation matrix
            if MATCHED_CC < 0.5:
                subcorrel += [-1]
            else:
            #subcorrel += [correlationCoefficient] --> correlation matrix
            # A block of code was removed for experiment *******************
                subcorrel += [BEST_MATCH]
        correlation += [subcorrel]
        print("Correlation: row " + str(c_row) + "/" + str(varianceHeight - 1))
    print("COMPLETED.")
    return correlation


def Residual(image, variance, window, varianceWidth, varianceHeight):
    """Function that calculate the residual of pixels in each window of the image\n
    Return a matrix of lists of residual of each pixel, in the size of variance matrix"""
    residual = []
    for v_row in range(varianceHeight):
        sub_resid = []
        for v_col in range(varianceWidth):
            # if it is not a point of interest, ignore the residual calculation
            if variance[v_row][v_col] == [0,0,0]:
                
                # check the surrounding points, if any of 8 nearby points is a distinct points in the original variance function *************
                # check this point as well. This is an experimental algorithm  --> edge expansion *****************
                if v_row == 0:
                    if v_col == 0:
                        # x x x
                        # x o
                        # x 
                        check = [variance[v_row][v_col+1], variance[v_row+1][v_col], variance[v_row+1][v_col+1]]
                    elif v_col == varianceWidth - 1:
                        # x x x
                        #   o x
                        #     x
                        check = [variance[v_row][v_col-1], variance[v_row+1][v_col], variance[v_row+1][v_col-1]]
                    else:
                        # x x x
                        #   o
                        # 
                        check = [variance[v_row][v_col-1], variance[v_row][v_col+1], variance[v_row+1][v_col], variance[v_row+1][v_col-1] , variance[v_row+1][v_col+1]]
                elif v_row == varianceHeight - 1:
                    if v_col == 0:
                        # x
                        # x o
                        # x x x
                        check = [variance[v_row][v_col+1], variance[v_row-1][v_col], variance[v_row-1][v_col+1]]
                    elif v_col == varianceWidth - 1:
                        #     x
                        #   o x
                        # x x x
                        check = [variance[v_row][v_col-1], variance[v_row-1][v_col], variance[v_row-1][v_col-1]]
                    else:
                        # 
                        #   o
                        # x x x
                        check = [variance[v_row][v_col-1], variance[v_row][v_col+1], variance[v_row-1][v_col], variance[v_row-1][v_col-1] , variance[v_row-1][v_col+1]]
                else:
                    if v_col == 0:
                        # x
                        # x o
                        # x
                        check = [variance[v_row][v_col+1], variance[v_row-1][v_col+1], variance[v_row-1][v_col], variance[v_row+1][v_col+1] , variance[v_row+1][v_col]]
                    elif v_col == varianceWidth - 1:
                        #     x
                        #   o x
                        #     x
                        check = [variance[v_row][v_col-1], variance[v_row-1][v_col-1], variance[v_row-1][v_col], variance[v_row+1][v_col-1] , variance[v_row+1][v_col]]
                    else:
                        # 
                        #   o
                        # 
                        check = [variance[v_row-1][v_col-1], variance[v_row-1][v_col], variance[v_row-1][v_col+1], variance[v_row][v_col-1], variance[v_row][v_col+1] , variance[v_row+1][v_col-1], variance[v_row+1][v_col], variance[v_row+1][v_col+1]]
                if check.count([0,0,0]) == len(check):
                    sub_resid += [[]]
                    continue
            
            # if it is a point of interest, collect all residuals of pixels in the windowed image
            i_left = v_col * window
            i_right = i_left + window
            i_top = v_row * window
            i_bottom = i_top + window
            # avg stores the average value of rgb in a window
            # rgb_window collects residuals of all pixels in a list
            r_rgb_avg = [0,0,0]
            r_rgb_window = []
            for i_row in range(i_top, i_bottom, 1):
                for i_col in range(i_left, i_right, 1):
                    r_rgb_subwindow = []
                    for c in range(3):
                        r_rgb_avg[c] += image[i_row][i_col][c]
                        r_rgb_subwindow += [image[i_row][i_col][c]]
                    r_rgb_window += [r_rgb_subwindow]
            # calculate averages of each channel in a window
            r_rgb_avg[0] = r_rgb_avg[0]  / (window**2)
            r_rgb_avg[1] = r_rgb_avg[1]  / (window**2)
            r_rgb_avg[2] = r_rgb_avg[2]  / (window**2)
            # calculate residuals of rgb channels for each pixel in a window
            for rgb in r_rgb_window:
                for c in range(3):
                    rgb[c] -= r_rgb_avg[c]
            sub_resid += [r_rgb_window]
        residual += [sub_resid]
        print("Residual: row " + str(v_row) + "/" + str(varianceHeight - 1))
    print("COMPLETED.")
    return residual


def Fill(corMatrix):
    """This function perform linear filling in both horizontal and vertical direction with linear fitting
    Return a correlation matrix that is filled"""
    filledMatrix = deepcopy(corMatrix)
    height = len(corMatrix)
    width = len(corMatrix[0])
    groundDispairty = 0
    for row in range(height):
        col = 0
        while (col < width):
            
            if row == height - 1:
                if corMatrix[row][col] > groundDispairty:
                    groundDispairty = corMatrix[row][col]
            
            headPixel = 0
            tailPixel = 0
            headDisparity = 0
            tailDisparity = 0
            if corMatrix[row][col] == -1:
                if col == 0:
                    filledMatrix[row][0] = 1
                    if corMatrix[row][1] != -1:
                        col += 1
                        continue
                    else:
                        headPixel = 0
                        headDisparity = 1
                        col = 1
                elif col == width - 1:
                    filledMatrix[row][col] = corMatrix[row][col-1]
                    break
                else:
                    headPixel = col - 1
                    headDisparity = corMatrix[row][headPixel]
                
                tailPixel = col + 1
                while (corMatrix[row][tailPixel] == -1):
                    tailPixel += 1
                    if tailPixel == width:
                        tailPixel = width - 1
                        filledMatrix[row][tailPixel] = 1
                        tailDisparity = 1
                        break
                if tailPixel != width - 1:
                    tailDisparity = corMatrix[row][tailPixel]
                
                for i in range(1, tailPixel-headPixel, 1):
                    filledMatrix[row][headPixel+i] = headDisparity + i * (tailDisparity-headDisparity)/(tailPixel-headPixel)
                col = tailPixel
            col += 1
        print("Fill: row " + str(row) + "/" + str(height - 1))
    print("COMPLETED.")
    
    for col in range(width):
        row = 0
        while (row < height):
            headPixel = 0
            tailPixel = 0
            headDisparity = 0
            tailDisparity = 0
            if corMatrix[row][col] == -1:
                if row == 0:
                    filledMatrix[0][col] = round((filledMatrix[0][col] + 1.1)/2)
                    if corMatrix[1][col] != -1:
                        row += 1
                        continue
                    else:
                        headPixel = 0
                        headDisparity = 1
                        row = 1
                elif row == height - 1:
                    filledMatrix[row][col] = round((filledMatrix[row][col] + groundDispairty)/2)
                    break
                else:
                    headPixel = row - 1
                    headDisparity = corMatrix[headPixel][col]
                
                tailPixel = row + 1
                while (corMatrix[tailPixel][col] == -1):
                    tailPixel += 1
                    if tailPixel == height:
                        tailPixel = height - 1
                        filledMatrix[tailPixel][col] = round((filledMatrix[tailPixel][col] + groundDispairty)/2)
                        tailDisparity = filledMatrix[tailPixel][col]
                        break
                if tailPixel != height - 1:
                    tailDisparity = corMatrix[tailPixel][col]
                
                for i in range(1, tailPixel-headPixel, 1):
                    filledMatrix[headPixel+i][col] += headDisparity + i * (tailDisparity-headDisparity)/(tailPixel-headPixel)
                    filledMatrix[headPixel+i][col] = round(filledMatrix[headPixel+i][col]/2)
                row = tailPixel
            row += 1
        print("Fill: col " + str(col) + "/" + str(width - 1))
    print("COMPLETED.")
    return filledMatrix


windowWidth = 5
fileData = ReadFile(0, windowWidth)

matplot.figure()
matplot.imshow(fileData[0])
matplot.figure()
matplot.imshow(fileData[1])

var_l = Variance(fileData[0], windowWidth)
matplot.figure()
matplot.imshow(var_l)
var_r = Variance(fileData[1], windowWidth)
matplot.figure()
matplot.imshow(var_r)

cor = Correlation(fileData[0], fileData[1], var_l, var_r, windowWidth, fileData[2])
matplot.figure()
matplot.imshow(cor)
filled_cor = Fill(cor)
matplot.figure()
matplot.imshow(filled_cor)

x = range(0, len(filled_cor), 1)
y = range(0, len(filled_cor[0]), 1)
xl = []
for x_e in x:
    for y_e in y:
        xl += [x_e]
yl = []
for x_e in x:
    for y_e in y:
        yl += [y_e]
zl = []
for x_e in x:
    for y_e in y:
        if (filled_cor[x_e][y_e] != 0):
            zl += [176.252 / filled_cor[x_e][y_e] * 209.059]
        else:
            zl += [176.252 * 209.059]

fig = matplot.figure()
ax = fig.add_subplot(111, projection='3d')
pnt3d = ax.scatter(xl,yl,zl, c=zl)
cbar = matplot.colorbar(pnt3d)
cbar.set_label("Distance (mm)")
matplot.show()