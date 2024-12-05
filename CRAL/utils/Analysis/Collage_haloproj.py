from PIL import Image, ImageEnhance

W  = 1920
H  = 640
 
halonum_arr = [1,2,4] 
for ihalo in halonum_arr :
  imgname_proj = 'proj_temp_noneq_z10.0_' + str(ihalo)+ '.png'
  imgname_tn   = 'tn_noneq_z10.0_' + str(ihalo)+ '.png'

  new = Image.new("RGBA", (int(W),int(H)))
  img1 = Image.open(imgname_proj)
  img2 = Image.open(imgname_tn)

  img1 = img1.resize((int(W/2),int(H)))
  img2 = img2.resize((int(W/2),int(H)))

  new.paste(img1, (0,0))
  new.paste(img2, (int(W/2),0))

  fignanme= 'halo_noneq_'+ str(ihalo) + '.png'
  print("saving :", fignanme)
  new.save(fignanme)