module MidiTS
    using DelimitedFiles, LinearAlgebra

    include("constr_series.jl")


    export get_time_series, get_pitch_seq
    export pitches_gaps, notas_hertz, get_onon_notes, serie_notas, get_new_voice, get_onoff_notes, silence_gaps, min_voces, max_tempo, rounding, filter_undef
    export indice, find_more_voices, get_subs

end
