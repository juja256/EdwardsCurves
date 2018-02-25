

def check_point_on_curve(cur, x, y):
    return (x**2 + y**2) % cur['p'] == (1 + cur['d'] * x**2 * y**2) % cur['p']


def main():
    cur192 = {
        'p': 0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEFFFFFFFFFFFFFFFF,
        'x': 0x44F083BB00E51AD91A2743284D31F57EE5C84826FCC91F4B,
        'y': 0x15FC16E5870524E0DBBE9EC8BB9F066C02A02B1978D4E029,
        'n': 0x3FFFFFFFFFFFFFFFFFFFFFFFEA75D4027230DD4DFFDB0455,
        'd': 0x6DBA6A,
    }

    print check_point_on_curve(cur192, 0x251CCE0526844F2963B9D3CB32EC6C3817AA7CED7BF3953B, 0xF8F3CE9E4CD4EF2866A3000DBED5CD4F2201332DC1F2D6E8)

if __name__ == '__main__':
    main()