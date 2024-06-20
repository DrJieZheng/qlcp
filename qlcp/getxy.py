#encoding=utf8

"""
click stars in the image, and get their x,y coordinates
"""


import matplotlib.pyplot as plt
import numpy as np
import astropy.io.fits as fits
from astropy.stats import sigma_clipped_stats
import pltgui
import datetime


def getxy(imgfile, size=10):
    print(datetime.datetime.now())
    # load image, and get clipped mean and std and image size
    img = fits.getdata(imgfile)
    imean, imed, istd = sigma_clipped_stats(img)
    ny, nx = img.shape
    
    # Construct the xy weight array for calculating barycentre
    wx = np.array([np.arange(nx)] * ny)
    wy = np.array([np.arange(ny)] * nx).T

    # Compute the barycentre according to x0,y0
    def barycentre(x0, y0, box=size, nit=5):
        # tranform to pixel coordinate
        x0, y0 = int(round(x0)), int(round(y0))
        # find the max value and revise the center
        for k in range(nit):
            # search box
            x0a, x0b = x0-box, x0+box
            y0a, y0b = y0-box, y0+box
            # find the max value and its index
            y1, x1 = np.unravel_index(np.argmax(img[y0a:y0b, x0a:x0b]), (box*2,box*2))
            x1, y1 = x0a+x1, y0a+y1
            # print(f"Max{k} | ({x0:4d},{y0:4d}) --> ({x1:7.2f},{y1:7.2f})  {img[y1,x1]:5d}")
            # 如果当前中心已经是最大值，就不迭代了
            # break if current center is the max value
            if x1 == x0 and y1 == y0:
                break
            # update center for next iteration
            x0, y0 = x1, y1
        
        # 检查最大值是否超过 5 sigma，如果没超过，说明点击的是背景
        # check if the max value is above 5 sigma, if not, it's a background
        if img[y0, x0] < imean + 5 * istd:
            x0, y0 = None, None
        else:
            for k in range(nit):
                # 分多轮计算重心，实际上每次只会跳1格，不知道为什么，但是足够了
                # search the barycentre in a box
                x0a, x0b = x0-box, x0+box
                y0a, y0b = y0-box, y0+box
                w = np.sum(img[y0a:y0b, x0a:x0b])
                x1 = np.sum(img[y0a:y0b, x0a:x0b] * wx[y0a:y0b, x0a:x0b]) / w
                y1 = np.sum(img[y0a:y0b, x0a:x0b] * wy[y0a:y0b, x0a:x0b]) / w
                x2, y2 = int(round(x1)), int(round(y1))
                # print(f"Cen{k} | ({x0:4d},{y0:4d}) --> ({x1:7.2f},{y1:7.2f})  {img[y2,x2]:5d}")
                # break if the barycentre is close to the center
                if x2 == x0 and y2 == y0:
                    break
                x0, y0 = x2, y2
        return x0, y0

    # 画图
    # plot the image
    fig = plt.figure(figsize=(8, 8))
    ax = fig.add_axes([0.05, 0.05, 0.9, 0.9])
    ax.imshow(img, vmin=imean-3*istd, vmax=imean+3*istd, cmap="gray", origin="lower")
    ax.set_title(imgfile)
    # 选中的目标标记，后续直接数组中增减
    # mark the selected targets, add/remove from array
    seltxt = [] # 每个元素是一个text对象 item is a text object
    selmar = [] # 每个元素是marker对象 item is a marker object
    selxy = []  # 每个元素是(x,y) item is (x,y)
    
    def point_sel(x, y):
        # 点击图像中之后找到点并标注
        # 先找重心，或者说就是现场找星，如果报告没找到，忽略
        # find the barycentre of star clicked, if not found, ignore it
        xx, yy = barycentre(x, y)
        if xx is None:
            return
        # 检查是否在其附近已经选了一个目标，如果是，就是删除，否则就是新增
        # compute 2D distance between selected points and the clicked point
        if len(selxy) > 0:
            dx = np.array([xx - xy[0] for xy in selxy])
            dy = np.array([yy - xy[1] for xy in selxy])
            dxy = np.sqrt(dx*dx+dy*dy)
            i = np.argmin(dxy)
        else:
            dxy = [size * 10]
            i = 0
        # if 2D distance is less than the given size, delete it
        if dxy[i] < size:
            # 删除
            seltxt[i].remove()
            del seltxt[i], selmar[i], selxy[i]
        else:
            # 新增
            seltxt.append(ax.text(xx+30, yy+30, "", 
                color="r", va="center", ha="center"))
            selmar.append(ax.scatter([xx], [yy], 
                color="none", marker="o", edgecolors="r"))
            selxy.append((xx, yy))
        # 重新标注标签文本
        # re-label the text
        for i, t in enumerate(seltxt):
            t.set_text(f"{i}")

    # 设置控制器
    # set the controller
    ctl = pltgui.i_btn_ctl(ax, image_action=point_sel)

    # 进入交互式
    # enter interactive mode
    plt.ion()
    plt.show()
    ctl.action_loop()
    plt.ioff()

    # 把结果进行格式化并输出
    # format and output the result
    xytxt = "\n".join(f"    ({xy[0]:4d},{xy[1]:4d})," for xy in selxy)
    print(f"{len(selxy)} points:")
    if len(selxy) > 0:
        print(f"[\n{xytxt}\n]")
    else:
        print("[]")
    
    return


if __name__ == "__main__":
    import zjarg
    txt, rval, rkey, other = zjarg.parse([])
    if txt:
        f = txt[0]
        size = int(rval[0]) if rval else 10
        getxy(f, size)
    else:
        print("请指定一个图像文件")
