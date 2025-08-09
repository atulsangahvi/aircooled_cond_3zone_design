# Streamlit application for air-cooled condenser design
import math
import streamlit as st
import pandas as pd
from CoolProp.CoolProp import PropsSI, get_global_param_string
from scipy.optimize import fsolve

st.set_page_config(layout="wide")
st.title("Air-Cooled Condenser Design Tool")

# ====================== INPUT SECTION ======================
with st.sidebar:
    st.header("🛠️ Geometry Inputs")
    tube_od_mm = st.number_input("Tube Outer Diameter (mm)", min_value=0.1, value=9.525)
    tube_thickness_mm = st.number_input("Tube Wall Thickness (mm)", min_value=0.01, value=0.35)
    row_pitch_mm = st.number_input("Row Pitch (mm)", min_value=0.1, value=25.4)
    tube_pitch_mm = st.number_input("Tube Pitch (mm)", min_value=0.1, value=25.4)
    fpi = st.number_input("Fins Per Inch", min_value=1, value=12)
    fin_thickness_mm = st.number_input("Fin Thickness (mm)", min_value=0.01, value=0.12)
    fin_material = st.selectbox("Fin Material", ["Aluminum", "Copper"])
    face_width_m = st.number_input("Coil Face Width (m)", min_value=0.1, value=1.0)
    face_height_m = st.number_input("Coil Face Height (m)", min_value=0.1, value=1.0)
    num_rows = st.number_input("Total Number of Rows", min_value=1, value=4)
    free_area_percent = st.slider("Free Flow Area %", 10, 100, 25)
    num_feeds = st.number_input("Number of Feeds", min_value=1, value=4)

    st.header("❄️ Refrigerant Inputs")
    fluid_list = get_global_param_string("FluidsList").split(',')
    refrigerants = sorted([f for f in fluid_list if f.startswith("R")])
    fluid = st.selectbox("Refrigerant", refrigerants, index=refrigerants.index("R134a") if "R134a" in refrigerants else 0)
    T1 = st.number_input("Inlet Superheat Temp (°C)", value=95.0)
    T3 = st.number_input("Outlet Subcooled Temp (°C)", value=52.0)
    T_cond = st.number_input("Condensing Temp (°C)", value=58.0)
    m_dot = st.number_input("Mass Flow Rate (kg/s)", min_value=0.01, value=0.6)
    air_temp = st.number_input("Air Inlet Temp (°C)", value=48.0)
    airflow_cmh = st.number_input("Air Flow (m³/hr)", min_value=1.0, value=10000.0)

# ====================== CALCULATIONS ======================
st.header("📊 Calculated Parameters")

# Derived Dimensions
tube_od_m = tube_od_mm / 1000
tube_thk_m = tube_thickness_mm / 1000
tube_id_m = max(0.001, tube_od_m - 2 * tube_thk_m)
tube_pitch_m = tube_pitch_mm / 1000
row_pitch_true_m = math.sqrt((tube_pitch_mm/2000)**2 + (row_pitch_mm/1000)**2)
fins_per_m = fpi * 39.37
fin_thk_m = fin_thickness_mm / 1000
fin_k = 235 if fin_material == "Aluminum" else 400
frontal_area_m2 = face_width_m * face_height_m
net_free_area = max(0.001, frontal_area_m2 * (free_area_percent/100))
airflow_m3s = max(0.001, airflow_cmh / 3600)
air_velocity_face = airflow_m3s / frontal_area_m2
air_velocity_fin = airflow_m3s / net_free_area

# Display geometry outputs
with st.expander("📐 Geometry Calculations", expanded=True):
    cols = st.columns(4)
    cols[0].metric("Tube OD", f"{tube_od_m*1000:.2f} mm")
    cols[1].metric("Tube ID", f"{tube_id_m*1000:.2f} mm")
    cols[2].metric("True Row Pitch", f"{row_pitch_true_m*1000:.1f} mm")
    cols[3].metric("Fins/m", f"{fins_per_m:.0f}")
    
    cols = st.columns(3)
    cols[0].metric("Frontal Area", f"{frontal_area_m2:.3f} m²")
    cols[1].metric("Free Flow Area", f"{net_free_area:.3f} m²")
    cols[2].metric("Face Velocity", f"{air_velocity_face:.2f} m/s")
    cols[0].metric("Fin Velocity", f"{air_velocity_fin:.2f} m/s")

# Tube and Fin Counts
tubes_per_row = max(1, math.floor(face_width_m / tube_pitch_m))
independent_circuits = max(1, math.floor(tubes_per_row / num_feeds))
m_dot_per_circuit = m_dot / independent_circuits
tube_length_per_tube = face_height_m
total_tubes = tubes_per_row * num_rows
total_tube_length = total_tubes * tube_length_per_tube
A_tube_ext = total_tube_length * math.pi * tube_od_m
A_tube_int = total_tube_length * math.pi * tube_id_m

# Fin Surface Area
num_fins = max(1, math.floor(face_height_m / 0.0254 * fpi))
L_fin = min(row_pitch_true_m, tube_pitch_m) / 2
A_fin_raw = row_pitch_true_m * num_rows * face_width_m * num_fins
A_hole = (math.pi/4 * tube_od_m**2) * total_tubes
A_fin = max(0, (A_fin_raw - A_hole) * 2)

# Air Properties
try:
    T_K = air_temp + 273.15
    mu_air = PropsSI("V", "T", T_K, "P", 101325, "Air")
    k_air = PropsSI("L", "T", T_K, "P", 101325, "Air")
    rho_air = PropsSI("D", "T", T_K, "P", 101325, "Air")
    cp_air = PropsSI("C", "T", T_K, "P", 101325, "Air")/1000  # kJ/kg-K
    Pr_air = PropsSI("PRANDTL", "T", T_K, "P", 101325, "Air")
    Re_air = max(1, rho_air*air_velocity_fin*tube_od_m/mu_air)
    C = 0.27 if Re_air < 1000 else 0.021
    n = 0.63 if Re_air < 1000 else 0.84
    Nu_air = C * Re_air**n * Pr_air**0.36
    h_air = Nu_air*k_air/tube_od_m if tube_od_m > 0 else 0
    C_air = max(0.001, airflow_m3s*rho_air*cp_air)  # kW/K
    
    with st.expander("🌬️ Air Properties", expanded=True):
        cols = st.columns(4)
        cols[0].metric("Density", f"{rho_air:.3f} kg/m³")
        cols[1].metric("Viscosity", f"{mu_air*1e6:.2f} μPa·s")
        cols[2].metric("Thermal Cond.", f"{k_air:.4f} W/m·K")
        cols[3].metric("Specific Heat", f"{cp_air:.3f} kJ/kg·K")
        
        cols = st.columns(4)
        cols[0].metric("Reynolds No.", f"{Re_air:.0f}")
        cols[1].metric("Prandtl No.", f"{Pr_air:.3f}")
        cols[2].metric("Nusselt No.", f"{Nu_air:.1f}")
        cols[3].metric("h (air side)", f"{h_air:.1f} W/m²K")
        
except Exception as e:
    st.error(f"❌ Error calculating air properties: {str(e)}")
    st.stop()

# Fin efficiency calculation
m = math.sqrt((2*h_air)/(fin_k*fin_thk_m)) if fin_k*fin_thk_m > 0 else 0
eta_fin = math.tanh(m*L_fin)/(m*L_fin) if m*L_fin > 0 else 0
A_air_total = A_tube_ext + A_fin*eta_fin
A_air_per_m = A_air_total/total_tube_length if total_tube_length > 0 else 0
total_tube_length_per_row = tubes_per_row * tube_length_per_tube

# Display surface areas
with st.expander("📏 Surface Areas", expanded=True):
    cols = st.columns(3)
    cols[0].metric("Tubes per Row", f"{tubes_per_row}")
    cols[1].metric("Independent Circuits", f"{independent_circuits}")
    cols[2].metric("Mass Flow per Circuit", f"{m_dot_per_circuit:.4f} kg/s")
    
    cols = st.columns(3)
    cols[0].metric("Total Tube Length", f"{total_tube_length:.1f} m")
    cols[1].metric("Tube Outer Area", f"{A_tube_ext:.2f} m²")
    cols[2].metric("Fin Area (effective)", f"{A_fin*eta_fin:.2f} m²")
    
    cols = st.columns(2)
    cols[0].metric("Total Air-Side Area", f"{A_air_total:.2f} m²")
    cols[1].metric("Area per Meter", f"{A_air_per_m:.2f} m²/m")

# Refrigerant Properties with zone-specific properties
try:
    P_cond = PropsSI("P", "T", T_cond + 273.15, "Q", 0, fluid)
    h1 = PropsSI("H", "P", P_cond, "T", T1 + 273.15, fluid)
    h2 = PropsSI("H", "P", P_cond, "Q", 1, fluid)
    h3 = PropsSI("H", "P", P_cond, "Q", 0, fluid)
    h4 = PropsSI("H", "P", P_cond, "T", T3 + 273.15, fluid)
    
    # Zone-specific average temperatures
    T_avg_desuper = (T1 + T_cond)/2
    T_avg_cond = T_cond
    T_avg_subcool = (T_cond + T3)/2
    
    # Zone-specific properties
    def get_ref_properties(P, T_avg, phase=None):
        if phase == "two-phase":
            return {
                'rho': PropsSI("D", "P", P, "Q", 0.5, fluid),
                'mu': PropsSI("V", "P", P, "Q", 0.5, fluid),
                'k': PropsSI("L", "P", P, "Q", 0.5, fluid),
                'cp': PropsSI("C", "P", P, "Q", 0.5, fluid)/1000,
                'Pr': PropsSI("PRANDTL", "P", P, "Q", 0.5, fluid)
            }
        else:
            return {
                'rho': PropsSI("D", "P", P, "T", T_avg + 273.15, fluid),
                'mu': PropsSI("V", "P", P, "T", T_avg + 273.15, fluid),
                'k': PropsSI("L", "P", P, "T", T_avg + 273.15, fluid),
                'cp': PropsSI("C", "P", P, "T", T_avg + 273.15, fluid)/1000,
                'Pr': PropsSI("PRANDTL", "P", P, "T", T_avg + 273.15, fluid)
            }
    
    props_desuper = get_ref_properties(P_cond, T_avg_desuper)
    props_cond = get_ref_properties(P_cond, T_avg_cond, "two-phase")
    props_subcool = get_ref_properties(P_cond, T_avg_subcool)
    
    with st.expander("🧊 Refrigerant Properties", expanded=True):
        cols = st.columns(4)
        cols[0].metric("Condensing Pressure", f"{P_cond/1000:.2f} kPa")
        cols[1].metric("Superheat Enthalpy", f"{h1/1000:.2f} kJ/kg")
        cols[2].metric("Vapor Enthalpy", f"{h2/1000:.2f} kJ/kg")
        cols[3].metric("Liquid Enthalpy", f"{h3/1000:.2f} kJ/kg")
        
        cols = st.columns(3)
        cols[0].metric("Subcool Enthalpy", f"{h4/1000:.2f} kJ/kg")
        cols[1].metric("Cp (desuperheat)", f"{props_desuper['cp']:.3f} kJ/kg·K")
        cols[2].metric("Cp (condensing)", f"{props_cond['cp']:.3f} kJ/kg·K")
        cols[0].metric("Cp (subcooling)", f"{props_subcool['cp']:.3f} kJ/kg·K")
        
except Exception as e:
    st.error(f"❌ Error calculating refrigerant properties: {str(e)}")
    st.stop()

# Heat Loads
Q_sens = m_dot * (h1-h2)/1000  # kW
Q_lat = m_dot * (h2-h3)/1000   # kW
Q_sub = m_dot * (h3-h4)/1000   # kW

# Calculate refrigerant-side heat transfer coefficients
def calculate_h_ref(props, m_dot_per_circuit, tube_id_m, T_avg, zone_type):
    """Calculate refrigerant-side h using appropriate correlations"""
    A_cross = (math.pi/4) * tube_id_m**2
    vel = m_dot_per_circuit / (props['rho'] * A_cross)
    Re = props['rho'] * vel * tube_id_m / props['mu']
    Pr = props['Pr']
    
    if zone_type == "condensing":
        # Use Akers et al. correlation for condensation
        G = m_dot_per_circuit / A_cross
        Re_eq = G * tube_id_m / props['mu']
        Nu = 0.026 * Re_eq**0.8 * Pr**0.33
    else:
        # For single-phase regions
        if Re < 2300:
            # Laminar flow - use Sieder-Tate correlation
            Nu = 1.86 * (Re * Pr * (tube_id_m/tube_length_per_tube))**0.33
        else:
            # Turbulent flow - Dittus-Boelter
            n = 0.4 if zone_type == "desuperheat" else 0.3  # heating or cooling
            Nu = 0.023 * Re**0.8 * Pr**n
    
    h_ref = Nu * props['k'] / tube_id_m
    return h_ref, Re, Pr, vel

h_ref_desuper, Re_desuper, Pr_desuper, vel_desuper = calculate_h_ref(
    props_desuper, m_dot_per_circuit, tube_id_m, T_avg_desuper, "desuperheat")
h_ref_cond, Re_cond, Pr_cond, vel_cond = calculate_h_ref(
    props_cond, m_dot_per_circuit, tube_id_m, T_avg_cond, "condensing")
h_ref_subcool, Re_subcool, Pr_subcool, vel_subcool = calculate_h_ref(
    props_subcool, m_dot_per_circuit, tube_id_m, T_avg_subcool, "subcooling")

# Calculate overall heat transfer coefficients
def calculate_U(h_ref, h_air, fin_efficiency, tube_k=400):
    """Calculate overall heat transfer coefficient with all resistances"""
    R_air = 1 / (h_air * fin_efficiency * A_air_total)
    R_ref = 1 / (h_ref * A_tube_int)
    R_wall = math.log(tube_od_m/tube_id_m)/(2*math.pi*tube_k*total_tube_length)
    R_fouling = 0.0001  # 10% fouling factor
    R_total = R_air + R_ref + R_wall + R_fouling
    return 1 / max(0.001, R_total)

U_desuper = calculate_U(h_ref_desuper, h_air, eta_fin)
U_cond = calculate_U(h_ref_cond, h_air, eta_fin)
U_subcool = calculate_U(h_ref_subcool, h_air, eta_fin)

with st.expander("🔥 Heat Transfer Parameters", expanded=True):
    cols = st.columns(4)
    cols[0].metric("Desuperheat h_ref", f"{h_ref_desuper:.1f} W/m²K")
    cols[1].metric("Condensing h_ref", f"{h_ref_cond:.1f} W/m²K")
    cols[2].metric("Subcooling h_ref", f"{h_ref_subcool:.1f} W/m²K")
    cols[3].metric("Air-side h", f"{h_air:.1f} W/m²K")
    
    cols = st.columns(3)
    cols[0].metric("Desuperheat U", f"{U_desuper:.1f} W/m²K")
    cols[1].metric("Condensing U", f"{U_cond:.1f} W/m²K")
    cols[2].metric("Subcooling U", f"{U_subcool:.1f} W/m²K")
    
    cols = st.columns(3)
    cols[0].metric("Desuperheat Q", f"{Q_sens:.2f} kW")
    cols[1].metric("Condensing Q", f"{Q_lat:.2f} kW")
    cols[2].metric("Subcooling Q", f"{Q_sub:.2f} kW")

# Flow regime diagnostics
with st.expander("🔍 Flow Regime Diagnostics", expanded=False):
    cols = st.columns(3)
    cols[0].metric("Desuperheat Re", f"{Re_desuper:.0f}", 
                  "Turbulent" if Re_desuper >= 2300 else "Laminar")
    cols[1].metric("Condensing Re_eq", f"{Re_cond:.0f}")
    cols[2].metric("Subcooling Re", f"{Re_subcool:.0f}", 
                  "Turbulent" if Re_subcool >= 2300 else "Laminar")
    
    cols = st.columns(3)
    cols[0].metric("Desuperheat Vel", f"{vel_desuper:.2f} m/s")
    cols[1].metric("Condensing Vel", f"{vel_cond:.2f} m/s")
    cols[2].metric("Subcooling Vel", f"{vel_subcool:.2f} m/s")

# ====================== ZONE CALCULATIONS ======================
def calculate_effectiveness(NTU, C_star):
    """Calculate effectiveness for crossflow with one fluid mixed"""
    if C_star == 0:  # Phase change
        return 1 - math.exp(-NTU)
    elif C_star == 1:
        return NTU / (1 + NTU)
    else:
        numerator = 1 - math.exp(-NTU*(1 - C_star))
        denominator = 1 - C_star*math.exp(-NTU*(1 - C_star))
        return numerator / denominator

zone_data = [
    ("Subcooling", Q_sub, U_subcool, h_ref_subcool, T3, props_subcool['cp']),
    ("Condensing", Q_lat, U_cond, h_ref_cond, T_cond, props_cond['cp']), 
    ("Desuperheat", Q_sens, U_desuper, h_ref_desuper, T1, props_desuper['cp'])
]

zone_outputs = []
diagnostics = []
air_t = air_temp
remaining_rows = num_rows
total_Q_actual = 0
total_Q_required = Q_sens + Q_lat + Q_sub

for label, Q_zone, U, h_ref, T_ref, cp_zone in zone_data:
    if Q_zone <= 0 or remaining_rows <= 0:
        zone_outputs.append((label, Q_zone, 0, 0, 0, air_t))
        diagnostics.append((label, 0, 0, 0, 0, 0, U, A_air_per_m, 0))
        continue
        
    # Capacity rates
    if label == "Condensing":
        C_min = C_air
        C_star = 0
        delta_T = max(1, T_ref - air_t)
    else:
        C_r = m_dot * cp_zone
        C_min = min(C_air, C_r)
        C_max = max(C_air, C_r)
        C_star = C_min / max(0.001, C_max)
        delta_T = max(1, abs(T_ref - air_t))

    # Calculate maximum possible heat transfer
    Q_max = C_min * delta_T
    
    # Numerically solve for NTU that gives required Q_zone
    NTU_guess = 1.0
    try:
        def equation(NTU):
            eps = calculate_effectiveness(NTU, C_star)
            return eps * Q_max - Q_zone
        
        NTU = max(0.1, min(fsolve(equation, NTU_guess)[0], 100))
        eps = calculate_effectiveness(NTU, C_star)
    except:
        NTU = 3.0  # Fallback value
        eps = min(0.99, Q_zone / Q_max)

    A_required = NTU * C_min * 1000 / max(1, U)
    tube_length_needed = A_required / max(0.1, A_air_per_m)
    rows_needed = tube_length_needed / max(0.1, total_tube_length_per_row)
    
    # Distribute rows proportionally to heat load
    rows_used = min(rows_needed, remaining_rows * (Q_zone/total_Q_required))
    rows_used = min(rows_used, remaining_rows)
    remaining_rows -= rows_used
    
    # Calculate actual heat transfer
    A_used = rows_used * total_tube_length_per_row * A_air_per_m
    Q_actual = eps * Q_max if rows_used >= rows_needed else U * A_used * delta_T / 1000
    total_Q_actual += Q_actual
    
    # Update air temperature (with physical limits)
    air_out = min(85, max(air_temp, air_t + (Q_actual/max(0.001, C_air))))
    
    zone_outputs.append((label, Q_zone, A_required, rows_needed, rows_used, air_out))
    diagnostics.append((label, C_min, C_star, delta_T, eps, NTU, U, A_air_per_m, Q_actual))
    
    air_t = air_out

# Final verification
if abs(total_Q_actual - total_Q_required) > 1:
    st.warning(f"Heat balance mismatch: {total_Q_actual:.2f} kW vs {total_Q_required:.2f} kW")

# ====================== OUTPUT SECTION ======================
st.header("📋 Design Results")

# Main summary
cols = st.columns(4)
cols[0].metric("Total Heat Load", f"{total_Q_required:.2f} kW")
cols[1].metric("Actual Heat Transfer", f"{total_Q_actual:.2f} kW")
cols[2].metric("Rows Utilization", f"{num_rows - remaining_rows:.2f}/{num_rows}")
cols[3].metric("Max Air Out Temp", f"{air_t:.1f} °C")

# Detailed results
st.subheader("🔍 Zone-by-Zone Results")
df = pd.DataFrame(zone_outputs, columns=[
    "Zone", "Q (kW)", "Area Needed (m²)", "Rows Needed", "Rows Used", "Air Out Temp (°C)"
])
st.dataframe(df.style.format({
    "Q (kW)": "{:.2f}",
    "Area Needed (m²)": "{:.2f}",
    "Rows Needed": "{:.2f}",
    "Rows Used": "{:.2f}",
    "Air Out Temp (°C)": "{:.1f}"
}))

# Diagnostic data
st.subheader("🛠️ Diagnostic Parameters")
diag_df = pd.DataFrame(diagnostics, columns=[
    "Zone", "C_min (kW/K)", "C*", "ΔT (°C)", "ε", "NTU", "U (W/m²K)", 
    "Area/m (m²/m)", "Q_actual (kW)"
])
st.dataframe(diag_df.style.format({
    "C_min (kW/K)": "{:.3f}",
    "C*": "{:.3f}",
    "ΔT (°C)": "{:.1f}",
    "ε": "{:.3f}",
    "NTU": "{:.2f}",
    "U (W/m²K)": "{:.1f}",
    "Area/m (m²/m)": "{:.2f}",
    "Q_actual (kW)": "{:.2f}"
}))

# Design verification
if (num_rows - remaining_rows) > num_rows * 1.01:
    st.error(f"⚠️ Requires {num_rows - remaining_rows:.2f} rows but only {num_rows} provided")
    st.warning(f"""
    **Recommendations:**
    - Increase number of rows (current: {num_rows})
    - Improve air-side heat transfer (current h_air: {h_air:.1f} W/m²K)
    - Check refrigerant flow rates (current: {m_dot:.3f} kg/s)
    - Verify temperature differences (min ΔT: {min([d[3] for d in diagnostics if d[3]>0]):.1f}°C)
    """)
elif sum([z[2] for z in zone_outputs]) > A_air_total * 1.05:
    st.warning(f"⚠️ Requires {sum([z[2] for z in zone_outputs]):.1f} m² but only {A_air_total:.1f} provided")
else:
    st.success("✅ Design meets requirements")

# Download buttons
csv1 = df.to_csv(index=False)
csv2 = diag_df.to_csv(index=False)

cols = st.columns(2)
cols[0].download_button(
    label="Download Zone Results",
    data=csv1,
    file_name="condenser_zone_results.csv",
    mime="text/csv"
)
cols[1].download_button(
    label="Download Diagnostics",
    data=csv2,
    file_name="condenser_diagnostics.csv",
    mime="text/csv"
)