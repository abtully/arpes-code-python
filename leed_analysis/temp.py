import os
import matplotlib.pyplot as plt

# ddir = r'E:\atully\LEED data\2022_February\28_Feb_2022'
# fn = r'enhanced_02_28_22_C60_graphene_SiC_1100am_24eV.tiff'
fname = os.path.abspath(r'E:\atully\LEED data\2022_February\28_Feb_2022\enhanced_02_28_22_C60_graphene_SiC_1100am_24eV.tiff')
# fp = os.path.join(ddir, fn)
# img_data = Image.open(fp)
print('about to load')
img_data = plt.imread(fname)
plt.imshow(img_data)
plt.show()

print('success')