"""Calculation example taken from EN 16282-1 Annex C.

Determination of the airflows for the kitchen of a canteen.
"""
from hvac import Quantity
from hvac.kitchen_ventilation import (
    KitchenAppliance,
    KitchenBlock,
    KitchenBlockArrangement,
    Kitchen,
    SupplyAirFlowType
)

Q_ = Quantity

# KITCHEN BLOCK 1: Kitchen appliances close together.

appliance1 = KitchenAppliance(
    name='electric cooker',
    P=Q_(16.9, 'kW'),
    q_dot_sen=Quantity(200, 'W / kW'),
    m_dot_w=Q_(118, 'g / (h * kW)'),
    h_d=Q_(1.2, 'm')
)

appliance2 = KitchenAppliance(
    name='electric deep fryer',
    P=Q_(14.4, 'kW'),
    q_dot_sen=Quantity(90, 'W / kW'),
    m_dot_w=Q_(1030, 'g / (h * kW)'),
    h_d=Q_(1.2, 'm')
)

appliance3 = KitchenAppliance(
    name='electric pressure cooker',
    P=Q_(12.4, 'kW'),
    q_dot_sen=Quantity(40, 'W / kW'),
    m_dot_w=Q_(15, 'g / (h * kW)'),
    h_d=Q_(1.2, 'm')
)

appliance4 = KitchenAppliance(
    name='electric tilting frying pan',
    P=Q_(14.4, 'kW'),
    q_dot_sen=Quantity(450, 'W / kW'),
    m_dot_w=Q_(588, 'g / (h * kW)'),
    h_d=Q_(1.2, 'm')
)

appliance5 = KitchenAppliance(
    name='electric pressure boiling pan 150 l',
    P=Q_(23.6, 'kW'),
    q_dot_sen=Quantity(40, 'W / kW'),
    m_dot_w=Q_(15, 'g / (h * kW)'),
    h_d=Q_(1.2, 'm')
)

appliance6 = KitchenAppliance(
    name='electric boiling pan 80 l',
    P=Q_(19.0, 'kW'),
    q_dot_sen=Quantity(35, 'W / kW'),
    m_dot_w=Q_(441, 'g / (h * kW)'),
    h_d=Q_(1.2, 'm')
)

appliance7 = KitchenAppliance(
    name='electric boiling pan 2x 60l',
    P=Q_(36.0, 'kW'),
    q_dot_sen=Quantity(35, 'W / kW'),
    m_dot_w=Q_(441, 'g / (h * kW)'),
    h_d=Q_(1.2, 'm')
)

appliance8 = KitchenAppliance(
    name='electric boiling pan 40 l',
    P=Q_(15.0, 'kW'),
    q_dot_sen=Quantity(35, 'W / kW'),
    m_dot_w=Q_(441, 'g / (h * kW)'),
    h_d=Q_(1.2, 'm')
)

kitchen_block1 = KitchenBlock(
    name='hood H1',
    width=Q_(1750, 'mm'),
    length=Q_(4250, 'mm'),
    arrangement=KitchenBlockArrangement.ANYWHERE,
    hooded=True
)
kitchen_block1.add_kitchen_appliances(
    appliance1,
    appliance2,
    appliance3,
    appliance4,
    appliance5,
    appliance6,
    appliance7,
    appliance8
)


# KITCHEN BLOCK 2: Kitchen appliances not close together.

appliance1 = KitchenAppliance(
    name='electric combination oven',
    P=Q_(36.5, 'kW'),
    q_dot_sen=Q_(120, 'W / kW'),
    m_dot_w=Q_(265, 'g / (h * kW)'),
    h_d=Q_(1.21, 'm'),
    width=Q_(895, 'mm'),
    length=Q_(885, 'mm')
)

appliance2 = KitchenAppliance(
    name='electric combination oven',
    P=Q_(9.3, 'kW'),
    q_dot_sen=Q_(120, 'W / kW'),
    m_dot_w=Q_(265, 'g / (h * kW)'),
    h_d=Q_(0.92, 'm'),
    width=Q_(895, 'mm'),
    length=Q_(885, 'mm')
)

appliance3 = KitchenAppliance(
    name='electric high pressure steamer',
    P=Q_(22.5, 'kW'),
    q_dot_sen=Q_(25, 'W / kW'),
    m_dot_w=Q_(294, 'g / (h * kW)'),
    h_d=Q_(1.3, 'm'),
    width=Q_(810, 'mm'),
    length=Q_(500, 'mm')
)

kitchen_block2 = KitchenBlock(
    name='hood H2',
    arrangement=KitchenBlockArrangement.AGAINST_WALL,
    hooded=True
)
kitchen_block2.add_kitchen_appliances(
    appliance1,
    appliance2,
    appliance3
)


# KITCHEN BLOCK 3: Unhooded kitchen appliances.

kitchen_block3 = KitchenBlock(
    name='unhooded',
    arrangement=KitchenBlockArrangement.AGAINST_WALL,
    hooded=False
)
kitchen_block3.add_kitchen_appliances(KitchenAppliance(
    name='electric refrigerator',
    P=Q_(0.5, 'kW'),
    q_dot_sen=Q_(700, 'W / kW'),
    m_dot_w=Q_(0.0, 'g / (h * kW)'),
    h_d=Q_(1.0, 'm'),
    width=Q_(800, 'mm'),
    length=Q_(800, 'mm'),
))


# KITCHEN

kitchen = Kitchen(
    a=SupplyAirFlowType.MIXED_FLOW_PLANE,
    f_simul=0.6
)
kitchen.add_kitchen_blocks(
    kitchen_block1,
    kitchen_block2,
    kitchen_block3
)
V_dot_ext_tot, V_dot_ext_unhooded, V_dot_ext_comp = kitchen.calc_extract_airflow()

print(
    f"total extract airflow kitchen = {V_dot_ext_tot:~P.0f}",
    f"extract air flow hood H1 = {kitchen_block1.V_dot_ext:~P.0f}",
    f"extract air flow hood H2 = {kitchen_block2.V_dot_ext:~P.0f}",
    f"extract air flow unhooded block= {kitchen_block3.V_dot_ext:~P.0f}",
    f"compensation airflow kitchen = {V_dot_ext_comp:~P.0f}",
    f"unhooded airflow kitchen = {V_dot_ext_unhooded:~P.0f}",
    sep='\n', end='\n\n'
)

print(
    "supply airflow directly blown into the hood of kitchen block 2 = "
    f"{kitchen_block2.V_dot_sup_dir:~P.3f}",
    sep='\n'
)
