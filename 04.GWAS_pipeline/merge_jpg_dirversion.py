#--coding:utf-8--
# merge jpg

from os import listdir
import sys
#pip3 install pillow
from PIL import Image

if __name__ == "__main__":
    
    nlist=[fn.split("_manhattan")[0] for fn in listdir(sys.argv[1]) if fn.endswith("manhattan.jpg")]

    
    for tn in nlist:
        im_list=[]
        im_list.append(Image.open(sys.argv[1]+"/"+tn+"_manhattan.jpg"))
        im_list.append(Image.open(sys.argv[1]+"/"+tn+"_QQ-Plot.jpg"))
        # 图片转化为相同的尺寸
        # ims = []
        # for i in im_list:
        #     #new_img = i.resize((1280, 1280), Image.BILINEAR)
        #     ims.append(i)

        # 单幅图像尺寸
        # width, height = ims[0].size
        width=3985 #2657+1328
        height=1328
        # 创建空白长图
        im_list[1]=im_list[1].resize((1328, 1328), Image.BILINEAR)
        result = Image.new(im_list[0].mode, (width, height ))
        

        # 拼接图片
        # for i, im in enumerate(ims):
        #     result.paste(im, box=(0, i * height))
        result.paste(im_list[0],box=(0,0))
        result.paste(im_list[1],box=(2657,0))

        # 保存图片
        result.save(sys.argv[1]+"/"+tn+"_merge.jpg")
