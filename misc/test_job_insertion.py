from ddglib.ppi_api import BindingAffinityDDGUserInterface

if __name__ == '__main__':
    ddGdb = BindingAffinityDDGUserInterface(passwd = read_file('ddgdb.pw').strip(), use_utf = False)
    ddG_interface = ddGInterface(passwd = read_file('ddgdb.pw').strip())

BindingAffinityDDGUserInterface(
