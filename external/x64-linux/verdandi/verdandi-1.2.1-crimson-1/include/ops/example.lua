-- Comment lines are introduced by "--".

-- An entry is in the form "name = value":
last_name = "Handel" -- for a string
birth_year = 1685    -- for a number
nationality = {"German", "English"} -- for a list of strings

-- Entries can be nested, as if there were placed in sections and subsections.
name = {
   first_name = "George",
   middle_name = "Frideric",
   last_name = "Handel"
}
compositions = {
   oratorios = {"Messiah", "Solomon", "Saul"},
   operas = {"Orlando", "Alcina", "Ariodante"},
   suites = {"Water Music", "Music for the Royal Fireworks"},
   concerti_grossi_op_6 = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12}
}

-- Thanks to Lua, it is possible to refer to any variable:
death_age = 1759 - birth_year
one_composition = compositions.suites[1] -- warning: indexes start at 1.
-- Concatenation of strings:
full_name = name.first_name.." "..name.last_name

-- Functions may be defined:
function sum(i, j, k)
   return i + j + k
end

function sum_product(i, j, k)
   return i + j + k, i * j * k
end
