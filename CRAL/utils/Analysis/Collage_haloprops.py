from PIL import Image, ImageEnhance

W  = 1920
H  = 640

new_sfr_vs_mstar   = Image.new("RGBA", (int(W),int(H)))

img1 = Image.open("sfr_vs_mstar_z10.png")
img2 = Image.open("sfr_vs_mstar_z9.png")
img3 = Image.open("sfr_vs_mstar_z8.png")

img1 = img1.resize((int(W/3),int(H)))
img2 = img2.resize((int(W/3),int(H)))
img3 = img3.resize((int(W/3),int(H)))

new_sfr_vs_mstar.paste(img1, (0,0))
new_sfr_vs_mstar.paste(img2, (int(W/3),0))
new_sfr_vs_mstar.paste(img3, (int(2*W/3),0))
figname = "collage_sfr_vs_mstar.png"
print("saving :",figname )
new_sfr_vs_mstar.save(figname)

'''
img4 = Image.open("sfr_vs_mhalo_z10.png")
img5 = Image.open("sfr_vs_mhalo_z9.png")
#img6 = Image.open("sfr_vs_mhalo_z10.png")
img4 = img4.resize((int(W/3),int(H)))
img5 = img5.resize((int(W/3),int(H)))
img6 = img6.resize((int(W/3),int(H)))

new_sfr_vs_mhalo.paste(img4, (0,0))
new_sfr_vs_mhalo.paste(img5, (int(W/3),0))
new_sfr_vs_mhalo.paste(img6, (int(2*W/3),0))
figname="collage_sfr_vs_mhalo.png"
print("saving :", figname)
new_sfr_vs_mhalo.save(figname)


img7 = Image.open("mstar_vs_mhalo_z12.png")
img8 = Image.open("mstar_vs_mhalo_z11.png")
img9 = Image.open("mstar_vs_mhalo_z10.png")

img7 = img7.resize((int(W/3),int(H)))
img8 = img8.resize((int(W/3),int(H)))
img9 = img9.resize((int(W/3),int(H)))

new_mstar_vs_mhalo.paste(img7, (0,0))
new_mstar_vs_mhalo.paste(img8, (int(W/3),0))
new_mstar_vs_mhalo.paste(img9, (int(2*W/3),0))
figname = "collage_mstar_vs_mhalo.png"
print("saving :", figname)
new_mstar_vs_mhalo.save(figname)


img10 = Image.open("mstar_to_mhalo_z12.png")
img11 = Image.open("mstar_to_mhalo_z11.png")
img12 = Image.open("mstar_to_mhalo_z10.png")

img10 = img10.resize((int(W/3),int(H)))
img11 = img11.resize((int(W/3),int(H)))
img12 = img12.resize((int(W/3),int(H)))

new_mstar_to_mhalo.paste(img10, (0,0))
new_mstar_to_mhalo.paste(img11, (int(W/3),0))
new_mstar_to_mhalo.paste(img12, (int(2*W/3),0))
figname = "collage_mstar_to_mhalo.png"
print("saving :", figname)
new_mstar_to_mhalo.save(figname)
'''