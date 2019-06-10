function pitches_gaps!(s, q) #Sometimes in the midifile there are small gaps between pitches, this function removes them.
    for v in s
        for i = 1:(size(v)[1]-1)
            m = abs(v[i,2] - v[i+1,1])
            #printlrintln(Voz[i,2],'\t',Voz[i+1,1],'\t',m)
            if m < q
                v[i,2] = v[i+1,1]
                #println(Voz[i,2],'\t', Voz[i+1,1])
            end
        end
      end
end
#########################################################################################################################################################################
function silence_gaps!(s,q) #Sometimes there are small gaps between silences and pitches, this function removes them.
    sm = q
    for v in s
        for i = 1:(size(v)[1]-1)
            m = mod(v[i,2] - v[i+1,1],sm)
            if m > 0
                v[i,2] = v[i,2] + (sm - m)
            end
        end
        x = length(v[:,1])
        n = mod(v[x,2] - v[x,1], sm)
        if n > 0
            v[x,2] = v[x,2] + (sm - n)
        end
    end
end
#################################################################################################################################################
function min_voces(Voces, div) #this function returns the minimum time value of all voices.
    ns = size(Voces)[1]
    mins = Array{Float64}(undef,ns)
    for i = 1:ns
        dif = filter(x->x >= div,Voces[i][:,2] - Voces[i][:,1])
        mins[i] = minimum(dif[find(dif)])
    end
    return minimum(mins)
end
################################################################################################################################################
function max_tempo(Voces, ns) #This function is to find what is the maximum time when the piece ends.
    Tfinal = Array{Float64}(undef,ns)
    for i = 1:ns
        Tfinal[i] = maximum(Voces[i][:,2])
    end
    return maximum(Tfinal)
end
################################################################################################################################
#rounding the possible fractional numbers obtained by the subdivision.
function rounding!(voces, q)
    nv = size(voces)[1]
    for i =1:nv
        voces[i][:,1] = map(x-> ceil(x/q), voces[i][:,1])
        voces[i][:,2] = map(x-> ceil(x/q), voces[i][:,2])
    end
end
####################################################################################################
##########################################################################################################
#Next functin filter the undef voices in left in the array of voices...
function filter_undef!(voces)
    nv = size(voces)[1]
    tt = Array{Bool}(undef,nv)
    for i=1:nv
        tt[i] = isassigned(voces,i)
    end
    deleteat!(voces,find(x->x==false,tt))
end
####################################################################################################################
function indice(tmax::Int64) #esta funcion solo regresa el arreglo que lleva el indice(numero) de las notas (o el eje x)
    ind = Array{Float64}(undef,tmax)
    for i = 1:tmax
        ind[i] = convert(Float64, i)
    end
    return(ind)
end
##############################################################################################################################
function notas_hertz!(notas::Array{Float64,1}) #esta funcion convierte el arreglo de entrada de numero de nota a su valor en hertz con A4 = 440Hz
    a4 = 440.0
    x = 1.0/12.0
    r = 2 ^ x
    for i = 1:length(notas)
        if notas[i] == 0.0 #si la nota ya es silencio, la omite
            continue
        end
        d = notas[i] - 69 # 69 es el A4
        notas[i] = a4 * r ^ d
    end
end
################################################################################
#Next function is to get the subgroups of size sz in a sequence, subgroups are with overlap.
function get_subs(t_s,sz)
    subs = []
    t = length(t_s)
    for i = 1:(t-sz+1)
        sub = []
        for j = i:(i+sz-1)
            push!(sub,t_s[j])
        end
        push!(subs,sub)
    end

    return subs
end
###################################################################################
#New functions to construct time series and pitch sequences
function get_onoff_notes(s)
    ini = s[s[:,3].==" Note_on_c", :] #getting the the starting places
    fin = s[s[:,3].==" Note_off_c",:] #ending
    dn = unique(ini[:,5]) #different notes.
    ndn = length(dn) #number of different notes
    posin_dn = Array{Vector}(undef, ndn) #initialize arrays for getting information of each note
    posfn_dn = Array{Vector}(undef, ndn)
    in_fi = Array{Matrix}(undef,ndn)
    for i = 1:ndn       #for each different note it gets when are they played
        tmp = ini[ini[:,5].==dn[i],:][:,2]
        tmp2 = fin[fin[:,5].==dn[i],:][:,2]
        l1 = length(tmp); l2 = length(tmp2)
        if l1 != l2
            nt = min(l1,l2)
            note = [dn[i] for j=1:nt]
            in_fi[i] = [tmp[1:nt] tmp2[1:nt] note]
        else
            note = [dn[i] for j = 1:length(tmp)]
            in_fi[i] = [tmp tmp2 note]
        end
    end
    nmat = vcat(in_fi...) #getting all notes in the same array
    return nmat[sortperm(nmat[:,1]),:]  #returns the array sorted by appereance.
end
################################################################################
function get_new_voice(v)
    new_v = []
    ids = []
    iids = [i for i=1:size(v)[1]]
    c = 1
    while c <= size(v)[1]-1
        if v[c+1,1] <= v[c,1]
            #println(v[c,:],'\t',v[c+1,:])
            push!(new_v, v[c+1,:])
            push!(ids,c+1)
            c+=1
        else
            #push!(ids,c)
            c+=1
        end
        #c += 1
        #println(c)
    end
    deleteat!(iids, ids)
    old_v = v[iids,:]
    new_v = convert(Array{Any,2},transpose(reshape(vcat(new_v...),3,:)))
    return [old_v, new_v]
end
function find_more_voices(v)
    v_1 =[]
    nvs = get_new_voice(v)
    push!(v_1,nvs[1])
    lv = nvs[2]
    br = false
    while br == false #this is to check if there are more voices
        ns = get_new_voice(lv)
        push!(v_1, ns[1])
        br = isempty(ns[2])
        if br; break; end
        #println(br)
        lv = ns[2]
    end
    return v_1
end
################################################################################
function get_onon_notes(s)
    ini = s[s[:,6].!=0, :] #getting the the starting places
    fin = s[s[:,6].==0,:] #ending
    dn = unique(ini[:,5]) #different notes.
    ndn = length(dn) #number of different notes
    posin_dn = Array{Vector}(undef, ndn) #initialize arrays for getting information of each note
    posfn_dn = Array{Vector}(undef, ndn)
    in_fi = Array{Matrix}(undef,ndn)
    for i = 1:ndn       #for each different note it gets when are they played
        tmp = ini[ini[:,5].==dn[i],:][:,2]
        tmp2 = fin[fin[:,5].==dn[i],:][:,2]
        l1 = length(tmp); l2 = length(tmp2)
        if l1 != l2
            nt = min(l1,l2)
            note = [dn[i] for j=1:nt]
            in_fi[i] = [tmp[1:nt] tmp2[1:nt] note]
        else
            note = [dn[i] for j = 1:length(tmp)]
            in_fi[i] = [tmp tmp2 note]
        end
    end
    nmat = vcat(in_fi...) #getting all notes in the same array
    return nmat[sortperm(nmat[:,1]),:]  #returns the array sorted by appereance.
end
################################################################################
#next function gets the multivariate time series of a pieces from a csv file
function get_time_series(s)
    sd = 16
    #mq / sd
    #sd is the smallest reference duration, it would divide the quarter of a note.
    mq = s[1,6] #quarter of a note
    nv = s[findlast(s[:,3], " Note_on_c"),1] - 1 #estimates how many voices are in the midi
    voces = Array{Matrix}(undef,nv) #initialze an array
    #next function is to get the time in miliseconds when the pitch starts and ends.
    if nv ==0
        nv = 1
        if findfirst(s[:,3], " Note_off_c") == 0
            b = s[s[:,1].==1,:] #takes the events of the voice i
            mat = b[b[:,3].==" Note_on_c", :]
            mat = mat[find(x->x!="", mat[:,5]),:]
            voces = get_onon_notes(mat) #construct an array of of information of initial time, finish time, pitch and intensity.
        else
            mat = s[s[:,1].==1,:]
            mat = mat[find(x->x!="", mat[:,5]),:]
            voces = get_onoff_notes(mat)
        end
    else
        voces = Array{Matrix}(undef,nv) #initialze an array
        if findfirst(s[:,3], " Note_off_c") == 0
            #checks if the midi has events of note_off
            for i = 2:(nv+1) #if does not, it construct the series in this way
                b = s[s[:,1].==i,:] #takes the events of the voice i
                if size(b[b[:,3].==" Note_on_c", :])[1] == 0; continue; end #checks if there are notes in the channel
                mat = b[b[:,3].==" Note_on_c", :]
                mat = mat[find(x->x!="", mat[:,5]),:]
                voces[i-1] =  get_onon_notes(mat)#construct an array of information of initial time, finish time, pitch and intensity.
            end
        else
            for i = 2:(nv+1) #if has note_off events , it does this way
                ini = findfirst(s[s[:,1].==i,3], " Note_on_c") #initial time
                fin = findlast(s[s[:,1].==i,3], " Note_off_c") #finish time
                if ini == 0 || fin == 0; continue; end
                mat = s[s[:,1].==i,:]
                mat = mat[find(x->x!="", mat[:,5]),:]
                voces[i-1] = get_onoff_notes(mat) #construct an array of information of initial time, finish time, pitch and intensity
            end
        end
    end
    if nv > 1
        filter_undef!(voces)
        filter!(x->length(x)>0,voces)
        n_vs = []
        nv = length(voces)
        for i=1:nv
            m_v = find_more_voices(voces[i])
            for j =1:length(m_v)
                push!(n_vs, m_v[j])
            end
        end
    else
        try
            n_vs = []
            m_v = find_more_voices(voces[1])
            for j =1:length(m_v)
                push!(n_vs, m_v[j])
            end
        catch
            n_vs = []
            m_v = find_more_voices(voces)
            for j =1:length(m_v)
                push!(n_vs, m_v[j])
            end
        end
    end
    filter!(x->length(x)>0,n_vs)
    nv = size(n_vs)[1] #the real number of voices
    q = mq/ round(Int, mq / min_voces(n_vs, mq / sd)) #defines the unit of time  q
    pitches_gaps!(n_vs,q) #fixes some gaps of miliseconds between the notes
    silence_gaps!(n_vs,q) #fixes some gaps of milisecons of rests
    rounding!(n_vs, q) #rounds up some possible decimals in the time
    tmax = round(Int,max_tempo(n_vs,nv)) #gets the total lenght of the piece
    series = zeros(tmax,nv+1)  #this is the initialization of the output array
    series[:,1] = indice(tmax) #this is just the index
    for i=2:(nv+1) #constructs the output time series
        series[:,i] = serie_notas(n_vs[i-1],tmax)
    end
    return series, q, mq
end
################################################################################
#get the sequences of pitches
function get_pitch_seq(s)
    sd = 16
    #mq / sd
    #sd is the smallest reference duration, it would divide the quarter of a note.
    mq = s[1,6] #quarter of a note
    nv = s[findlast(s[:,3], " Note_on_c"),1] - 1 #estimates how many voices are in the midi
    voces = Array{Matrix}(undef,nv) #initialze an array
    #next function is to get the time in miliseconds when the pitch starts and ends.
    if nv ==0
        nv = 1
        if findfirst(s[:,3], " Note_off_c") == 0
            b = s[s[:,1].==1,:] #takes the events of the voice i
            mat = b[b[:,3].==" Note_on_c", :]
            mat = mat[find(x->x!="", mat[:,5]),:]
            voces = get_onon_notes(mat) #construct an array of of information of initial time, finish time, pitch and intensity.
        else
            mat = s[s[:,1].==1,:]
            mat = mat[find(x->x!="", mat[:,5]),:]
            voces = get_onoff_notes(mat)
        end
    else
        voces = Array{Matrix}(undef,nv) #initialze an array
        if findfirst(s[:,3], " Note_off_c") == 0
            #checks if the midi has events of note_off
            for i = 2:(nv+1) #if does not, it construct the series in this way
                b = s[s[:,1].==i,:] #takes the events of the voice i
                if size(b[b[:,3].==" Note_on_c", :])[1] == 0; continue; end #checks if there are notes in the channel
                mat = b[b[:,3].==" Note_on_c", :]
                mat = mat[find(x->x!="", mat[:,5]),:]
                voces[i-1] =  get_onon_notes(mat)#construct an array of information of initial time, finish time, pitch and intensity.
            end
        else
            for i = 2:(nv+1) #if has note_off events , it does this way
                ini = findfirst(s[s[:,1].==i,3], " Note_on_c") #initial time
                fin = findlast(s[s[:,1].==i,3], " Note_off_c") #finish time
                if ini == 0 || fin == 0; continue; end
                mat = s[s[:,1].==i,:]
                mat = mat[find(x->x!="", mat[:,5]),:]
                voces[i-1] = get_onoff_notes(mat) #construct an array of information of initial time, finish time, pitch and intensity
            end
        end
    end
    if nv > 1
        filter_undef!(voces)
        filter!(x->length(x)>0,voces)
        n_vs = []
        nv = length(voces)
        for i=1:nv
            m_v = find_more_voices(voces[i])
            for j =1:length(m_v)
                push!(n_vs, m_v[j])
            end
        end
    else
        try
            n_vs = []
            m_v = find_more_voices(voces[1])
            for j =1:length(m_v)
                push!(n_vs, m_v[j])
            end
        catch
            n_vs = []
            m_v = find_more_voices(voces)
            for j =1:length(m_v)
                push!(n_vs, m_v[j])
            end
        end
    end
    filter!(x->length(x)>0,n_vs)
    nv = size(n_vs)[1] #the real number of voices
    pitch_seq = Array{Vector,1}(undef, nv)
    for i = 1:nv
        pitch_seq[i] = n_vs[i][:,3]
    end
    return pitch_seq
end
