pro sav2json,sav_fname
; assumes IDL standard name.sav
names=STRSPLIT(sav_fname, '.', ESCAPE='\', /EXTRACT)
json_fname=names[0]+'.json'
restore,sav_fname,/v
sim=JSON_SERIALIZE(sims,PRECISION=3,/pretty)
openw,1,json_fname
printf,1,sim
close,1
print,'all data are now saves into:',json_fname
end
