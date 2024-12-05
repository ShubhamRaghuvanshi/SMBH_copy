from PIL import Image, ImageEnhance

W  = 1680
H  = 1740
 
new = Image.new("RGBA", (int(W),int(H)))

img1 = Image.open("gasdens_proj_star+sink_z12.png")
img2 = Image.open("gasdens_proj_star+sink_z11.png")
img3 = Image.open("gasdens_proj_star+sink_z10.png")

img4 = Image.open("gasdens_proj_star+sn+sink_z12.png")
img5 = Image.open("gasdens_proj_star+sn+sink_z11.png")
img6 = Image.open("gasdens_proj_star+sn+sink_z10.png")

img7 = Image.open("gasdens_proj_star+sink+agn_z12.png")
img8 = Image.open("gasdens_proj_star+sink+agn_z11.png")
img9 = Image.open("gasdens_proj_star+sink+agn_z10.png")

img10 = Image.open("gasdens_proj_star+sn+sink+agn_z12.png")
img11 = Image.open("gasdens_proj_star+sn+sink+agn_z11.png")
img12 = Image.open("gasdens_proj_star+sn+sink+agn_z10.png")


img1 = img1.resize((int(W/3),int(H/4)))
img2 = img2.resize((int(W/3),int(H/4)))
img3 = img3.resize((int(W/3),int(H/4)))

img4 = img4.resize((int(W/3),int(H/4)))
img5 = img5.resize((int(W/3),int(H/4)))
img6 = img6.resize((int(W/3),int(H/4)))

img7 = img7.resize((int(W/3),int(H/4)))
img8 = img8.resize((int(W/3),int(H/4)))
img9 = img9.resize((int(W/3),int(H/4)))

img10 = img10.resize((int(W/3),int(H/4)))
img11 = img11.resize((int(W/3),int(H/4)))
img12 = img12.resize((int(W/3),int(H/4)))


new.paste(img1, (0,0))
new.paste(img2, (int(W/3),0))
new.paste(img3, (int(2*W/3),0))

new.paste(img4, (0,int(H/4)))
new.paste(img5, (int(W/3),int(H/4)))
new.paste(img6, (int(2*W/3),int(H/4)))

new.paste(img7, (0,int(2*H/4)))
new.paste(img8, (int(W/3),int(2*H/4)))
new.paste(img9, (int(2*W/3),int(2*H/4)))

new.paste(img10, (0,int(3*H/4)))
new.paste(img11, (int(W/3),int(3*H/4)))
new.paste(img12, (int(2*W/3),int(3*H/4)))

figname = "collage_density_proj.png"
print("saving :", figname)
new.save(figname)