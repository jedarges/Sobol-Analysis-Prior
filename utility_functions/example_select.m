function data = example_select(spec)

switch spec
    case 'line_fitting'
         data = line_fitting_data;
    case 'seir'
         data = seir_data;
end