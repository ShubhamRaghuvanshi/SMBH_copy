from PIL import Image, ImageEnhance

W  = 1920
H  = 640

new_msink_vs_mhalo          = Image.new("RGBA", (int(W),int(H))) 
new_msinktomstar_vs_mhalo   = Image.new("RGBA", (int(W),int(H))) 


img1 = Image.open("msink_vs_mhalo_z12.png")
img2 = Image.open("msink_vs_mhalo_z11.png")
img3 = Image.open("msink_vs_mhalo_z10.png")
img1 = img1.resize((int(W/3),int(H)))
img2 = img2.resize((int(W/3),int(H)))
img3 = img3.resize((int(W/3),int(H)))

new_msink_vs_mhalo.paste(img1, (0,0))
new_msink_vs_mhalo.paste(img2, (int(W/3),0))
new_msink_vs_mhalo.paste(img3, (int(2*W/3),0))
figname = "collage_msink_vs_mhalo.png"
print("saving :", figname)
new_msink_vs_mhalo.save(figname)

img4 = Image.open("msinktomstar_vs_mhalo_z12.png")
img5 = Image.open("msinktomstar_vs_mhalo_z11.png")
img6 = Image.open("msinktomstar_vs_mhalo_z10.png")
img4 = img4.resize((int(W/3),int(H)))
img5 = img5.resize((int(W/3),int(H)))
img6 = img6.resize((int(W/3),int(H)))

new_msinktomstar_vs_mhalo.paste(img4, (0,0))
new_msinktomstar_vs_mhalo.paste(img5, (int(W/3),0))
new_msinktomstar_vs_mhalo.paste(img6, (int(2*W/3),0))
figname="collage_msinktomstar_vs_mhalo.png"
print("saving :", figname)
new_msinktomstar_vs_mhalo.save(figname)
