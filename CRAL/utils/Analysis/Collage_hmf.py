from PIL import Image, ImageEnhance

W  = 1920
H  = 640
 
new = Image.new("RGBA", (int(W),int(H)))

img1 = Image.open("hmf_cmp_z12.png")
img2 = Image.open("hmf_cmp_z11.png")
img3 = Image.open("hmf_cmp_z10.png")


img1 = img1.resize((int(W/3),int(H)))
img2 = img2.resize((int(W/3),int(H)))
img3 = img3.resize((int(W/3),int(H)))


new.paste(img1, (0,0))
new.paste(img2, (int(W/3),0))
new.paste(img3, (int(2*W/3),0))

figname = "collage_hmf.png"
print("saving :", figname) 
new.save(figname)