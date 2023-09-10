#--coding:utf-8--
# merge jpg

from os import listdir
import sys
#pip3 install pillow
from PIL import Image

if __name__ == "__main__":
    
    nlist=[fn.split("_ldheatmap")[0] for fn in listdir(sys.argv[1]) if fn.endswith("_ldheatmap.jpg")]

    
    for tn in nlist:
        im_list=[]
        im_list.append(Image.open(sys.argv[1]+"/"+tn+"_manhattan.jpg"))
        im_list.append(Image.open(sys.argv[1]+"/"+tn+"_ldheatmap.jpg"))
        # 图片转化为相同的尺寸
        # ims = []
        # for i in im_list:
        #     #new_img = i.resize((1280, 1280), Image.BILINEAR)
        #     ims.append(i)

        # 单幅图像尺寸
        # width, height = ims[0].size
        width=2657+ 47+65#2567+1328
        height=1328+1805
        # 创建空白长图
        #im_list[1]=im_list[1].resize((1328, 1328), Image.BILINEAR)
        im_list[1] = im_list[1].crop((0, 900, 2704, 2704))
        result = Image.new(im_list[0].mode, (width, height ),(255,255,255))
        

        # 拼接图片
        # for i, im in enumerate(ims):
        #     result.paste(im, box=(0, i * height))
        result.paste(im_list[0],box=(0,0))
        result.paste(im_list[1],box=(65,1328))

        # 保存图片
        result.save(sys.argv[1]+"/"+tn+"_merge.jpg")
